# export HDF5_DISABLE_VERSION_CHECK=1
from ast import Constant, Expression
import numpy as np
from dolfin import *
from fenics import *
import time
import ufl

# -----------------------------------------------------------------------------
# Cauchy Stress Tensor
# -----------------------------------------------------------------------------
def sigma(u, p, vg_l, vg_c, vg_i, k, f0, s0, prm_HO):
    """
    Computes the Cauchy and Piola-Kirchhoff stress tensors for the uterus model.
    
    Parameters
    ----------
        u : Displacement field.
        p : Pressure term.
        vg_l, vg_c, vg_i : Growth parameters (longitudinal, circumferential, isotropic).
        k : Region indicator (0 = uterine body, 1 = cervix).
        f0, s0 : Local fiber directions (longitudinal and circumferential).
        prm_HO : Material parameters [c_cvx, c, c_1, c_2, zeta1].
    
    Returns
    -------
        Fe : Elastic deformation gradient.
        sigma_elast : Elastic Cauchy stress.
        sigma_tot : Total Cauchy stress (including active effects).
        C : Right Cauchy-Green deformation tensor.
    """

    d = len(u)
    I = Identity(d)
    F = I + grad(u)
    C = F.T * F
    J = det(F)

    # Growth tensor
    Fg = (k * I + (1.0 - k)* (vg_i * I 
                              + (vg_l - 1.0) * outer(f0, f0) 
                              + (vg_c - 1.0) * outer(s0, s0)))

    Fe = F * inv(Fg)
    Ce = Fe.T * Fe
    Je = det(Fe)

    # Fiber invariants
    I4 = inner(Ce * f0, f0)  # Longitudinal fibers
    I6 = inner(Ce * s0, s0)  # Circumferential fibers

    # Material parameters
    c_cvx, c, c_1, c_2, zeta1 = prm_HO

    # Subplus function (x+ = max(x, 0))
    subplus = lambda var: conditional(ge(var, 0.0), var, 0.0)

    # Second Piola-Kirchhoff stress tensor for elastic deformation
    Sel = (k * (c_cvx * I - p * Je * inv(Ce).T) 
        + (1.0 - k) * ( c * I - p * Je * inv(Ce).T 
                   + 2.0 * zeta1 * c_1 * subplus(I4 - 1.0) 
                   * exp(c_2 * (subplus(I4 - 1.0) ** 2)) * outer(f0, f0) 
                   + 2.0 * (1.0 - zeta1) * c_1 * subplus(I6 - 1.0) 
                   * exp(c_2 * (subplus(I6 - 1.0) ** 2)) * outer(s0, s0)))
    

    # Push-forward (contravariant)
    tauKirchE = Fe * Sel * Fe.T

    # Elastic contribution (Cauchy stress)
    sigma_elast = (1.0 / Je) * tauKirchE

    # Second Piola-Kirchhoff stress tensor
    S = inv(Fg) * Sel * inv(Fg).T

    # Push-forward (contravariant)
    tauKirch = F * S * F.T

    # Total Cauchy stress
    sigma_tot = (1.0 / J) * tauKirch

    return Fe, sigma_elast, sigma_tot, C

# -----------------------------------------------------------------------------
# Function to compute endometrial (uterine cavity) volume
# -----------------------------------------------------------------------------
def compute_intrauterine_volume(J, F, mesh, ds, u=None):
    """
    Computes the intrauterine (endometrial cavity) volume in liters (L).

    Parameters
    ----------
        J : Determinant of the deformation gradient (Jacobian).
        F : Deformation gradient tensor.
        mesh : Finite element mesh of the uterus model.
        ds : Boundary measure identifying labeled surfaces.
        u : Displacement field. If provided, the deformed configuration is used.

    Returns
    -------
        vol_L : Uterine cavity volume in liters (L).
    """

    X = SpatialCoordinate(mesh)
    n_mesh = FacetNormal(mesh)

    # Label 20 = Endometrium region
    if u is not None:
        vol_mm3 = abs(
            assemble((-1.0 / 3.0) * dot(X + u, J * inv(F).T * n_mesh) * ds(20))
        )  # Volume in mm3 (deformed configuration)
    else:
        vol_mm3 = abs(
            assemble((-1.0 / 3.0) * dot(X, n_mesh) * ds(20))
        )  # Volume in mm3 (reference configuration)

    # Convert mm3 -> liters (1 mm3 = 1e-6 mL = 1e-6e-3 L = 1e-9 L)
    vol_L = 1e-6 * vol_mm3

    return vol_L

# -----------------------------------------------------------------------------
# Function to compute mesh (uterine) volume
# -----------------------------------------------------------------------------
def compute_uterine_volume(u, dx):
    """
    Computes the volume of the deformed uterine tissue.

    Parameters:
    -----------
    displacement : Displacement field of the mesh.
    dx_measure : dolfin.Measure
        Integration measure (dx) over the mesh.

    Returns:
    --------
    float
        Volume of the uterus in liters (L).
    """
    # Calculate the deformation gradient
    deformation_grad = Identity(3) + grad(u)
    
    # Compute the volume by integrating the determinant of the deformation gradient
    vol_mm3 = assemble(det(deformation_grad)*dx)

    # Convert mm3 -> liters (1 mm3 = 1e-6 mL = 1e-6e-3 L = 1e-9 L)
    vol_L = 1e-6 * vol_mm3

    # Convert from mm3 to m3
    return vol_L


# -----------------------------------------------------------------------------
# Function to solves the uterine growth and mechanics problem using a
# variational formulation
# -----------------------------------------------------------------------------
def UterusSolve(fileOutput, file_results, mesh, bdry, materials, f0, s0, tinicial, 
                prm_Holzapfel, prm_GrowthL, prm_GrowthC, IUP, k_mola, alpha, nsteps):
    """
    The function simulates mechanical deformation and growth of the uterus under
    intrauterine pressure, using a nonlinear finite element approach with a mixed
    formulation (displacement + pressure + growth factors).

    Parameters
    ----------
    fileOutput : File for writing intermediate fields.
    file_results : File for writing solution fields at each time step.
    mesh : Computational mesh of the uterus.
    bdry : Boundary markers for Dirichlet/Neumann conditions.
    materials : Material markers for different regions.
    f0 : Longitudinal fiber direction.
    s0 : Circumferential fiber direction.
    tinicial : Initial simulation time.
    prm_Holzapfel : Material parameters for the Holzapfel model.
    prm_GrowthL : Parameters for longitudinal growth.
    prm_GrowthC : Parameters for circumferential growth.
    IUP : Intrauterine pressure.
    k_mola : Spring constants for different regions.
    alpha : Isotropic growth rate.
    nsteps : Number of time steps.

    Returns
    -------
    vol_intra : Intrauterine volume at final step.
    vol_utero : Uterine tissue volume at final step.
    """

    # -----------------------------------------------------------------------------
    # Variational formulation and solution
    # -----------------------------------------------------------------------------
    # Vector element - u:
    P2v = VectorElement('CG', mesh.ufl_cell(), 2)
    # Scalar element - p, vgl, vgc, vgi:
    P1 = FiniteElement('CG', mesh.ufl_cell(), 1)
    # Function space formed by elements (u,p,vg_l,vg_c,vg_i)
    W = FunctionSpace(mesh, MixedElement(P2v, P1, P1, P1, P1))

    # Create function space
    w = Function(W)
    dW = TrialFunction(W)
    psiW = TestFunction(W)

    # Split functions into displacement and other variables
    (u, p, vg_l, vg_c, vg_i) = split(w)
    (v, q, wg_l, wg_c, wg_i) = split(psiW)

    # Print mesh info
    print(f'\n\nNumber cells =  {mesh.num_cells()}')
    print(f'Number edges =  {mesh.num_edges()}')
    print(f'Number vertices = {mesh.num_vertices()}')
    print(f'Number dof = {W.tabulate_dof_coordinates().shape[0]}\n\n')
    
    # Define class for determining k
    class KExpression(UserExpression):
        def eval(self, values, x):
            if 90 <= x[0] <= 140.0:
                values[0] = 1.0
            else:
                values[0] = 0.0

        def value_shape(self):
            return ()
             
    # Define function space for k
    V0 = FunctionSpace(mesh, 'DG', 0)   
    
    # Create and interpolate function k
    k = KExpression(degree = 1)
    k_function = interpolate(k, V0)
    k_function.rename('k','k')
    fileOutput.write(k_function,0)    
        
    # Define new measures for interior domains and exterior boundaries
    dx = Measure('dx', domain=mesh, subdomain_data=materials)
    ds = Measure('ds',domain=mesh, subdomain_data=bdry)

    # Kinematics
    d = u.geometric_dimension()
    I = Identity(d)             # Identity tensor
    F = variable(I + grad(u))   # Deformation gradient
    J = det(F)

    # Growth tensor - Body: k=0, Cervix: k=1  
    Fg = (k*I + (1.0 - k)*(vg_i*I 
                          + (vg_l - 1.0) * outer(f0,f0) 
                          + (vg_c - 1.0) * outer(s0,s0)))
    
    # Elastic tensor Fe
    Fe = F*inv(Fg)
    # Elastic volume change
    Je = det(Fe)
    # Right Cauchy-Green tensor for elastic deformation
    Ce = Fe.T*Fe
    # Total right Cauchy-Green tensor (elastic + growth)
    C = F.T * F

    # Ce invariants
    I4 = inner(Ce*f0,f0)     # Longitudinal fiber invariant
    I6 = inner(Ce*s0,s0)     # Circumferential fiber invariant

    # Subplus function: returns variable if >=0, else 0
    subplus = lambda var : conditional(ge(var, 0.0), var, 0.0)
    
    # Material constants
    c_cvx, c, c_1, c_2, zeta1 = prm_Holzapfel 

    # Second Piola-Kirchhoff stress tensor for elastic deformation
    # Body: k=0, Cervix: k=1
    Se = (k * (c_cvx * I - p * Je * inv(Ce).T) 
        + (1.0 - k) * ( c * I - p * Je * inv(Ce).T 
                   + 2.0 * zeta1 * c_1 * subplus(I4 - 1.0) 
                   * exp(c_2 * (subplus(I4 - 1.0) ** 2)) * outer(f0, f0) 
                   + 2.0 * (1.0 - zeta1) * c_1 * subplus(I6 - 1.0) 
                   * exp(c_2 * (subplus(I6 - 1.0) ** 2)) * outer(s0, s0)))
    
    # Second Piola-Kirchhoff stress tensor
    S = inv(Fg)*Se*(inv(Fg).T)
    # First Piola-Kirchhoff stress tensor
    P = F*S
    # Incompressibility condition
    hydpress = Je - 1.0  
    
    # Boundary markers: peri-body=30, endo-body=20, base-body=10, peri-cervix=60, endo-cervix=40, base-cervix=50
    # Define new markers
    new_marker_post = 70
    new_marker_pelv = 80    
    new_marker_base = 90
    
    # Update facet markers based on position
    for facet in facets(mesh):
        x_face = facet.midpoint().x()
        y_face = facet.midpoint().y()
        rotulo_atual = bdry[facet]
        if (y_face < -5.0) and (rotulo_atual == 30 or rotulo_atual == 60):  # posterior
            bdry[facet] = new_marker_post
        elif (x_face > 127.0 and x_face < 132.0) and (rotulo_atual == 50 or rotulo_atual == 40):  # base
            bdry[facet] = new_marker_base
        elif (x_face > 80.0 and x_face < 140.0) and (y_face >= -5.0) and (rotulo_atual in [30, 50, 60]):  # inferior pelvic
            bdry[facet] = new_marker_pelv
    
    # Rename function for correct legend display
    bdry.rename("Label", "")          # ParaView will now show "Label" instead of "f"
    
    # Save modified facet markers
    File(f"facet_markers_modified.pvd") << bdry

    # Define Dirichlet boundary condition
    zerov = Constant((0.0,0.0,0.0))
    bc1 = DirichletBC(W.sub(0), zerov, bdry, 90)   # Fix base-cervix - marker 90
    bcs = [bc1]
    
    # Spring constants
    k_peri, k_peri_post, k_peri_pelv = k_mola
    
    # Step-wise loading (for plotting and convergence)
    tfinal = 33.57142857 - tinicial
    dt = Constant(tfinal/nsteps)
    pressure = np.linspace(0, IUP, nsteps)
 
    # -------------------------------------------------------------------------
    # Weak form definition for equilibrium and growth problem
    # -------------------------------------------------------------------------

    # Normal vector to mesh facets
    N = FacetNormal(mesh)

    # Initial applied pressure
    p0 = Constant(pressure[0])

    # Weak form of mechanical equilibrium
    eq1 = ( inner(P, grad(v))*dx               # Elastic stress contribution
            + dot(p0*J*inv(F).T*N, v)*ds(20)  # Endometrial pressure on label 20
            + dot(k_peri*u, v)*ds(30)         # Pericervical spring forces
            + dot(k_peri_post*u, v)*ds(70)    # Posterior spring forces
            + dot(k_peri_pelv*u, v)*ds(80)    # Pelvic spring forces
            + hydpress*q*dx )                 # Incompressibility constraint


    # -------------------------------------------------------------------------
    # Growth problem: initialize growth variables
    # -------------------------------------------------------------------------

    # Initial growth factors = 1
    expr0 = Expression('1.0', degree=3) 
    vg0_l = interpolate(expr0, W.sub(2).collapse())
    vg0_c = interpolate(expr0, W.sub(3).collapse())
    vg0_i = interpolate(expr0, W.sub(4).collapse())

    # Set initial vectors to 1
    w.sub(2).vector().set_local(np.ones(w.vector().size()))
    w.sub(3).vector().set_local(np.ones(w.vector().size()))
    w.sub(4).vector().set_local(np.ones(w.vector().size()))
    w.vector().apply("")

    # Growth parameters
    tau_l, vmax_l, lcrit_l, gamma_l = prm_GrowthL
    tau_c, vmax_c, lcrit_c, gamma_c = prm_GrowthC 

    # Iteration constants
    t = Constant(0.0)
    step = 0

    # Compute total and elastic stretches along fibers
    l_long = sqrt(inner(C*f0,f0))   # Total longitudinal stretch
    lg_l = vg_l       
    le_l = l_long/lg_l       
    l_circ = sqrt(inner(C*s0,s0))   # Total circumferential stretch
    lg_c = vg_c     
    le_c = l_circ/lg_c       

    # Heaviside-like function to smoothly activate growth
    def heav_tanh(x, lcrit):
        b = 150.0
        y = (1.0/2.0)*(1.0+(exp((2.0*b)*(x-lcrit))-1.0)/(exp((2.0*b)*(x-lcrit))+1.0))
        return y

    # Growth evolution equations
    eq_long = ((vg_l - vg0_l)*wg_l*dx 
               - dt*(1.0/tau_l) * (((vmax_l - vg_l)/(vmax_l - 1.0))**gamma_l)
               * (heav_tanh(le_l,lcrit_l))*wg_l*dx )
    eq_circ = ((vg_c - vg0_c)*wg_c*dx 
               - dt*(1.0/tau_c) * (((vmax_c - vg_c)/(vmax_c - 1.0))**gamma_c)
               * (heav_tanh(le_c,lcrit_c))*wg_c*dx )
    eq_iso = (vg_i - vg0_i)*wg_i*dx - dt*alpha*wg_i*dx

    # Total weak form including equilibrium and growth
    eq = eq1 + eq_long + eq_circ + eq_iso

    # Jacobian for Newton solver
    Jac = derivative(eq, w, dW)

    # -------------------------------------------------------------------------
    # Create and configure the nonlinear variational problem
    # -------------------------------------------------------------------------
    # Define the variational problem R(w) = 0
    problem = NonlinearVariationalProblem(eq, w, bcs, J=Jac)

    # Create a solver for the nonlinear problem
    solver  = NonlinearVariationalSolver(problem)

    # Solver parameters
    solver.parameters['nonlinear_solver']                    = 'newton'   # Use Newton's method
    solver.parameters['newton_solver']['linear_solver']      = 'mumps'    # Linear solver
    solver.parameters['newton_solver']['absolute_tolerance'] = 1e-5      # Absolute convergence
    solver.parameters['newton_solver']['relative_tolerance'] = 1e-5      # Relative convergence
    solver.parameters['newton_solver']['maximum_iterations'] = 50        # Max iterations

    # Record start time
    tempo_inicial = time.time() # in seconds
    
    # Clear output file before simulation
    with open(f'./results_measures/out_p13_t_press_volIU_volUter.txt', 'w') as f_out:
        f_out.write("")  

    # -------------------------------------------------------------------------
    # Time-stepping loop
    # -------------------------------------------------------------------------
    while (step < nsteps):
        tempo_inicial1 = time.time()
        print(f"\n--- SIMULATION ---")
        print(f"Step: {step}")

        # Update intrauterine pressure
        p0.assign(pressure[step])

        # Solve nonlinear variational problem for current step
        solver.solve()
        
        # Extract solution components (displacement, pressure, growth factors)
        uaux, paux, vglaux, vgcaux, vgiaux = w.split(deepcopy = True)
        uaux.rename("u","displacement")
        paux.rename("p","pressure")
        vglaux.rename('vg_long','growth-factor_long')
        vgcaux.rename('vg_circ','growth-factor_circ')
        vgiaux.rename('vg_iso','growth-factor_isot')

        # Update previous growth factors for next step
        vg0_l.assign(vglaux)
        vg0_c.assign(vgcaux)
        vg0_i.assign(vgiaux)

        # Write current step results to files
        file_results.write(uaux, step)
        file_results.write(paux, step)
        file_results.write(vglaux, step)
        file_results.write(vgcaux, step)
        file_results.write(vgiaux, step)

        # -------------------------------------------------------------------------
        # Compute stresses, fiber stretches, and other post-processing quantities
        # -------------------------------------------------------------------------

        # Compute elastic deformation gradient, elastic and total Cauchy stress, total C
        Fe, sig_e, sig_t, Ct = Sigma(uaux, paux, vglaux, vgcaux, vgiaux, k, f0, s0, prm_Holzapfel)  

        # Project and save elastic Cauchy stress
        sigmaE_proj = project(sig_e, TensorFunctionSpace(mesh, "DG", 0))
        sigmaE_proj.rename("Sigma", "Elastic Cauchy Stress Tensor")
        file_results.write(sigmaE_proj, step)

        # Project and save total Cauchy stress
        sigmaT_proj = project(sig_t, TensorFunctionSpace(mesh, "DG", 0))
        sigmaT_proj.rename("Sigma_t", "Total Cauchy Stress Tensor")
        file_results.write(sigmaT_proj, step)

        # Project and save elastic deformation gradient
        Fe_proj = project(Fe, TensorFunctionSpace(mesh, "DG", 0))
        Fe_proj.rename("Fe", "Elastic Deformation Gradient")
        file_results.write(Fe_proj, step)

        # -------------------------------------------------------------------------
        # Compute total fiber stretches along longitudinal and circumferential directions
        # -------------------------------------------------------------------------
        lambd_l = sqrt(inner(Ct*f0, f0))   # Longitudinal fiber stretch
        lambd_c = sqrt(inner(Ct*s0, s0))   # Circumferential fiber stretch

        l_l = project(lambd_l, FunctionSpace(mesh, "DG", 0))
        l_l.rename("lambda_l", "Longitudinal Fiber Stretch")
        print("lambda_long MIN MAX", l_l.vector().get_local().min(), l_l.vector().get_local().max())

        l_c = project(lambd_c, FunctionSpace(mesh, "DG", 0))
        l_c.rename("lambda_c", "Circumferential Fiber Stretch")
        print("lambda_circ MIN MAX", l_c.vector().get_local().min(), l_c.vector().get_local().max())

        file_results.write(l_l, step)
        file_results.write(l_c, step)

        # -------------------------------------------------------------------------
        # Compute unit vectors and fiber stresses
        # -------------------------------------------------------------------------
        F = I + grad(u)
        f_def = F*f0
        f_uni = f_def / sqrt(inner(f_def, f_def))   # Unit vector along deformed longitudinal fiber

        s_def = F*s0
        s_uni = s_def / sqrt(inner(s_def, s_def))   # Unit vector along deformed circumferential fiber

        # Stress along fibers (unit vectors)
        stress_fib_l = inner(sig_t*f_uni, f_uni)
        stress_fib_c = inner(sig_t*s_uni, s_uni)

        s_f_l = project(stress_fib_l, FunctionSpace(mesh, "DG", 0))
        s_f_l.rename("Stress_Fiber_l", "Stress Fiber Longitudinal")
        s_f_c = project(stress_fib_c, FunctionSpace(mesh, "DG", 0))
        s_f_c.rename("Stress_Fiber_c", "Stress Fiber Circumferential")

        file_results.write(s_f_l, step)
        file_results.write(s_f_c, step)

        # -------------------------------------------------------------------------
        # Compute Von-Mises stress (deviatoric)
        # -------------------------------------------------------------------------
        s_dev = sig_t - (1.0/3.0)*(tr(sig_t))*I
        v_mises = sqrt((3.0/2.0)*inner(s_dev, s_dev))

        V_mises_proj = project(v_mises, FunctionSpace(mesh, "DG", 0))
        V_mises_proj.rename("Von_Mises", "Von-Mises Stress")
        file_results.write(V_mises_proj, step)

        # -------------------------------------------------------------------------
        # Compute intrauterine and total uterine volumes, convert pressure to kPa
        # -------------------------------------------------------------------------
        pres_kpa = 1000.0 * p0
        vol_intra = compute_intrauterine_volume(J, F, mesh, ds, u)
        vol_utero = compute_uterine_volume(u, dx)

        tempo_final1 = time.time()
        print(f'Time: {float(t)} ')
        print(f'Volume intrauterino (L): {vol_intra} ')
        print(f'Volume total utero (L): {vol_utero} ')
        print(f'Pressure (kPa): {float(pres_kpa)} ')

        # Update time and step
        t.assign(t + dt)
        step += 1

        # Append current results to output file
        with open(f'./results_measures/out_p13_t_press_volIU_volUter.txt', 'a') as f_out:
            f_out.write(f"{float(t):.6e}\t{float(pres_kpa):.6e}\t{float(vol_intra):.10e}\t{float(vol_utero):.10e}\n")

        # Print iteration and total simulation time
        print(f"Iteration time: {(tempo_final1 - tempo_inicial1)/60.0} minutes")
        print(f"Total elapsed time: {(tempo_final1 - tempo_inicial)/60.0} minutes")

    # -------------------------------------------------------------------------
    # Close results file and print total execution time
    # -------------------------------------------------------------------------
    file_results.close()
    tempo_final = time.time()
    print(f"Total execution time: {(tempo_final - tempo_inicial)/60.0} minutes")