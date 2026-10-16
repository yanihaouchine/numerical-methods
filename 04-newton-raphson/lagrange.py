import numpy as np

def elastic_force(k, x0, y0):
    """
    Creates an elastic force (with zero natural length) centered at (x0, y0) with intensity k.
    Returns a function that computes the force and its Jacobian at a given point.
    """
    def force(x, y):
        # Force vector: f = k * ([x, y] - [x0, y0])
        fx = k * (x - x0)
        fy = k * (y - y0)
        
        # Jacobian matrix:
        # [∂fx/∂x ∂fx/∂y]
        # [∂fy/∂x ∂fy/∂y]
        dfx_dx = k
        dfx_dy = 0
        dfy_dx = 0
        dfy_dy = k
        
        return np.array([fx, fy]), np.array([[dfx_dx, dfx_dy], [dfy_dx, dfy_dy]])
    
    return force

def centrifugal_force(k, x0, y0):
    """
    Creates a centrifugal force centered at (x0, y0) with intensity k.
    Returns a function that computes the force and its Jacobian at a given point.
    """
    def force(x, y):
        # Force vector: f = k * ([x, y] - [x0, y0])
        fx = k * (x - x0)
        fy = k * (y - y0)
        
        # Jacobian matrix:
        # [∂fx/∂x ∂fx/∂y]
        # [∂fy/∂x ∂fy/∂y]
        dfx_dx = k
        dfx_dy = 0
        dfy_dx = 0
        dfy_dy = k
        
        return np.array([fx, fy]), np.array([[dfx_dx, dfx_dy], [dfy_dx, dfy_dy]])
    
    return force

def gravitational_force(k, x0, y0):
    """
    Creates a gravitational force centered at (x0, y0) with intensity k.
    Returns a function that computes the force and its Jacobian at a given point.
    """
    def force(x, y):
        dx = x - x0
        dy = y - y0
        r_squared = dx**2 + dy**2
        r_cubed = r_squared * np.sqrt(r_squared)
        
        # Guard against division by zero
        if r_cubed < 1e-10:
            r_cubed = 1e-10
        
        # Force vector: f = -k * ([x, y] - [x0, y0]) / |[x, y] - [x0, y0]|^3
        fx = -k * dx / r_cubed
        fy = -k * dy / r_cubed
        
        # Jacobian matrix calculations
        common_term = k / r_cubed
        r_fifth = r_cubed * r_squared
        
        # ∂fx/∂x = -k * (1/r^3 - 3*dx^2/r^5)
        dfx_dx = -common_term * (1 - 3 * dx**2 / r_squared)
        
        # ∂fx/∂y = k * 3*dx*dy/r^5
        dfx_dy = common_term * 3 * dx * dy / r_squared
        
        # ∂fy/∂x = k * 3*dx*dy/r^5
        dfy_dx = common_term * 3 * dx * dy / r_squared
        
        # ∂fy/∂y = -k * (1/r^3 - 3*dy^2/r^5)
        dfy_dy = -common_term * (1 - 3 * dy**2 / r_squared)
        
        return np.array([fx, fy]), np.array([[dfx_dx, dfx_dy], [dfy_dx, dfy_dy]])
    
    return force

def combine_forces(forces):
    """
    Combines multiple force functions into a single force function.
    
    Parameters:
    forces -- List of force functions that each return a force vector and Jacobian
    
    Returns:
    A function that computes the combined force and Jacobian
    """
    def total_force(x, y):
        total_f = np.zeros(2)
        total_jacobian = np.zeros((2, 2))
        
        for force_func in forces:
            f, jacobian = force_func(x, y)
            total_f += f
            total_jacobian += jacobian
            
        return total_f, total_jacobian
    
    return total_force

from scipy.linalg import solve

def newton_raphson(force_func, initial_guess, max_iterations=100, tolerance=1e-10):
    """
    Implements the Newton-Raphson method for finding equilibrium points.
    
    Parameters:
    force_func -- Function that returns force and its Jacobian
    initial_guess -- Initial position [x, y]
    max_iterations -- Maximum number of iterations
    tolerance -- Convergence criterion
    
    Returns:
    Equilibrium position, force at equilibrium, and convergence status
    """
    position = np.array(initial_guess, dtype=float)
    
    for i in range(max_iterations):
        force, jacobian = force_func(position[0], position[1])
        
        # Check if we're already at an equilibrium
        force_magnitude = np.linalg.norm(force)
        if force_magnitude < tolerance:
            return position, force, True, i
        
        # Newton-Raphson update: p_{n+1} = p_n - J^(-1) * f
        # We need to solve J * Δp = -f for Δp
        try:
            delta = solve(jacobian, -force)
            position += delta
        except np.linalg.LinAlgError:
            # Jacobian is singular, can't proceed
            return position, force, False, i
        
        # Check if the correction is small enough
        if np.linalg.norm(delta) < tolerance:
            # Get final force and Jacobian at the converged position
            force, jacobian = force_func(position[0], position[1])
            return position, force, True, i
    
    # Get final force and Jacobian at the non-converged position
    force, jacobian = force_func(position[0], position[1])
    return position, force, False, max_iterations

# Set up the specific forces for the problem
def setup_system():
    # Gravitational force with k=1 at (0,0)
    grav1 = gravitational_force(1.0, 0.0, 0.0)
    
    # Gravitational force with k=0.01 at (1,0)
    grav2 = gravitational_force(0.01, 1.0, 0.0)
    
    # Calculate the barycenter of the two masses
    m1 = 1.0  # Mass 1 (inferred from gravitational coefficient)
    m2 = 0.01  # Mass 2 (inferred from gravitational coefficient)
    total_mass = m1 + m2
    barycenter_x = (m1 * 0.0 + m2 * 1.0) / total_mass
    barycenter_y = (m1 * 0.0 + m2 * 0.0) / total_mass
    
    # Centrifugal force with k=1 at the barycenter
    centri = centrifugal_force(1.0, barycenter_x, barycenter_y)
    
    # Combine all forces
    system_force = combine_forces([grav1, grav2, centri])
    
    return system_force, barycenter_x, barycenter_y

# Verify the implementation with the test point
def verify_system():
    system_force, bx, by = setup_system()
    test_point = [1.5, 0.0]
    
    force, jacobian = system_force(test_point[0], test_point[1])
    
    print(f"Barycenter: ({bx}, {by})")
    print(f"At point {test_point}:")
    print(f"Force: {force}")
    print(f"Jacobian:\n{jacobian}")
    
    # Reference values
    ref_force = np.array([1.00565457, 0.0])
    ref_jacobian = np.array([[1.75259259, 0.0], [0.0, 0.6237037]])
    
    # Check if our calculations match the reference values
    force_error = np.linalg.norm(force - ref_force)
    jacobian_error = np.linalg.norm(jacobian - ref_jacobian)
    
    print(f"Force error: {force_error}")
    print(f"Jacobian error: {jacobian_error}")
    
    return force_error < 1e-6 and jacobian_error < 1e-6

# Find equilibrium points
def find_equilibrium_points():
    system_force, bx, by = setup_system()
    
    # Define initial guesses for different equilibrium points
    # These are educated guesses based on known Lagrangian points in similar systems
    initial_guesses = [
        [bx, 0.0],           # L1 - between the two masses
        [1.2, 0.0],          # Beyond the smaller mass
        [-0.5, 0.0],         # Beyond the larger mass
        [0.5, 0.5],          # Above the line connecting the masses
        [0.5, -0.5]          # Below the line connecting the masses
    ]
    
    equilibrium_points = []
    
    for i, guess in enumerate(initial_guesses):
        position, force, converged, iterations = newton_raphson(system_force, guess)
        
        if converged:
            equilibrium_points.append({
                'position': position,
                'force_magnitude': np.linalg.norm(force),
                'iterations': iterations,
                'initial_guess': guess
            })
            
            print(f"Equilibrium point #{len(equilibrium_points)}")
            print(f"  Position: ({position[0]:.8f}, {position[1]:.8f})")
            print(f"  Force magnitude: {np.linalg.norm(force):.8e}")
            print(f"  Iterations: {iterations}")
            print(f"  Initial guess: {guess}")
            
            # Get eigenvalues of the Jacobian to check stability
            _, jacobian = system_force(position[0], position[1])
            eigenvalues = np.linalg.eigvals(jacobian)
            print(f"  Eigenvalues: {eigenvalues}")
            print()
    
    return equilibrium_points
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Import your existing functions
# (I'm assuming these functions are already defined as in your code)

def plot_force_field():
    # Set up the system as in your code
    system_force, barycenter_x, barycenter_y = setup_system()
    
    # Define the grid for our vector field
    x = np.linspace(-0.5, 1.5, 40)
    y = np.linspace(-1.0, 1.0, 40)
    X, Y = np.meshgrid(x, y)
    
    # Calculate forces at each grid point
    U = np.zeros_like(X)
    V = np.zeros_like(Y)
    magnitude = np.zeros_like(X)
    
    for i in range(len(x)):
        for j in range(len(y)):
            force, _ = system_force(X[j, i], Y[j, i])
            mag = np.linalg.norm(force)
            
            # Normalize for better visualization
            if mag > 1e-10:
                U[j, i] = force[0] / mag
                V[j, i] = force[1] / mag
            else:
                U[j, i] = 0
                V[j, i] = 0
                
            magnitude[j, i] = mag
    
    # Find equilibrium points
    initial_guesses = [
        [barycenter_x, 0.0],  # L1 - between the two masses
        [1.2, 0.0],           # Beyond the smaller mass
        [-0.5, 0.0],          # Beyond the larger mass
        [0.5, 0.5],           # Above the line connecting the masses
        [0.5, -0.5]           # Below the line connecting the masses
    ]
    
    equilibrium_points = []
    for guess in initial_guesses:
        position, force, converged, _ = newton_raphson(system_force, guess)
        if converged and np.linalg.norm(force) < 1e-6:
            # Check if this point is already in our list (within some tolerance)
            is_duplicate = False
            for p in equilibrium_points:
                if np.linalg.norm(p['position'] - position) < 1e-4:
                    is_duplicate = True
                    break
            
            if not is_duplicate:
                # Get eigenvalues to determine stability
                _, jacobian = system_force(position[0], position[1])
                eigenvalues = np.linalg.eigvals(jacobian)
                
                # Classify stability
                if all(np.real(eigenvalues) > 0):
                    stability = "Unstable"
                    marker = 'X'
                    color = 'red'
                elif all(np.real(eigenvalues) < 0):
                    stability = "Stable"
                    marker = 'o'
                    color = 'green'
                else:
                    stability = "Saddle"
                    marker = 's'
                    color = 'orange'
                
                equilibrium_points.append({
                    'position': position,
                    'eigenvalues': eigenvalues,
                    'stability': stability,
                    'marker': marker,
                    'color': color
                })
    
    # Create figure and axis
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))
    
    # Plot 1: Force field vectors
    ax1.set_title('Force Field with Equilibrium Points')
    ax1.quiver(X, Y, U, V, magnitude, cmap='viridis', pivot='mid', alpha=0.8)
    
    # Plot masses
    mass1 = Circle((0, 0), 0.05, color='blue', label='Mass 1 (k=1.0)')
    mass2 = Circle((1, 0), 0.02, color='green', label='Mass 2 (k=0.01)')
    bary = Circle((barycenter_x, barycenter_y), 0.01, color='black', label='Barycenter')
    
    ax1.add_patch(mass1)
    ax1.add_patch(mass2)
    ax1.add_patch(bary)
    
    # Plot equilibrium points
    for i, eq in enumerate(equilibrium_points):
        pos = eq['position']
        ax1.plot(pos[0], pos[1], marker=eq['marker'], markersize=10, 
                 color=eq['color'], label=f"L{i+1}: {eq['stability']}")
    
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.axis('equal')
    ax1.grid(True)
    ax1.legend()
    
    # Plot 2: Force magnitude contour
    contour = ax2.contourf(X, Y, np.log10(magnitude + 1e-10), 50, cmap='viridis')
    ax2.set_title('Log Force Magnitude')
    
    # Add colorbar
    divider = make_axes_locatable(ax2)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(contour, cax=cax, label='Log10(Force Magnitude)')
    
    # Plot masses and equilibrium points again
    ax2.add_patch(Circle((0, 0), 0.05, color='blue'))
    ax2.add_patch(Circle((1, 0), 0.02, color='green'))
    ax2.add_patch(Circle((barycenter_x, barycenter_y), 0.01, color='black'))
    
    for eq in equilibrium_points:
        pos = eq['position']
        ax2.plot(pos[0], pos[1], marker=eq['marker'], markersize=10, color=eq['color'])
    
    ax2.set_xlabel('x')
    ax2.set_ylabel('y')
    ax2.axis('equal')
    ax2.grid(True)
    
    plt.tight_layout()
    
    # Print equilibrium point information
    print("Equilibrium Points:")
    for i, eq in enumerate(equilibrium_points):
        pos = eq['position']
        print(f"L{i+1} ({pos[0]:.6f}, {pos[1]:.6f}) - {eq['stability']}")
        print(f"  Eigenvalues: {eq['eigenvalues']}")
    
    return fig, equilibrium_points

# Function to plot the potential function (negative of the work done by forces)
def plot_potential():
    system_force, barycenter_x, barycenter_y = setup_system()
    
    # Define the grid for the potential calculation
    x = np.linspace(-0.5, 1.5, 200)
    y = np.linspace(-1.0, 1.0, 200)
    X, Y = np.meshgrid(x, y)
    
    # Calculate potential (we'll approximate it by integrating from infinity)
    potential = np.zeros_like(X)
    
    # Calculate gravitational potential
    # For mass 1 at (0,0)
    r1 = np.sqrt(X**2 + Y**2)
    potential -= 1.0 / r1  # k=1.0
    
    # For mass 2 at (1,0)
    r2 = np.sqrt((X-1)**2 + Y**2)
    potential -= 0.01 / r2  # k=0.01
    
    # Centrifugal potential (negative because it's repulsive)
    rc = np.sqrt((X-barycenter_x)**2 + (Y-barycenter_y)**2)
    potential += 0.5 * 1.0 * rc**2  # k=1.0
    
    # Find equilibrium points
    initial_guesses = [
        [barycenter_x, 0.0],  # L1 - between the two masses
        [1.2, 0.0],           # Beyond the smaller mass
        [-0.5, 0.0],          # Beyond the larger mass
        [0.5, 0.5],           # Above the line connecting the masses
        [0.5, -0.5]           # Below the line connecting the masses
    ]
    
    equilibrium_points = []
    for guess in initial_guesses:
        position, force, converged, _ = newton_raphson(system_force, guess)
        if converged and np.linalg.norm(force) < 1e-6:
            # Check if this point is already in our list
            is_duplicate = False
            for p in equilibrium_points:
                if np.linalg.norm(p['position'] - position) < 1e-4:
                    is_duplicate = True
                    break
            
            if not is_duplicate:
                # Get eigenvalues to determine stability
                _, jacobian = system_force(position[0], position[1])
                eigenvalues = np.linalg.eigvals(jacobian)
                
                # Classify stability
                if all(np.real(eigenvalues) > 0):
                    stability = "Unstable"
                    marker = 'X'
                    color = 'red'
                elif all(np.real(eigenvalues) < 0):
                    stability = "Stable"
                    marker = 'o'
                    color = 'green'
                else:
                    stability = "Saddle"
                    marker = 's'
                    color = 'orange'
                
                equilibrium_points.append({
                    'position': position,
                    'eigenvalues': eigenvalues,
                    'stability': stability,
                    'marker': marker,
                    'color': color
                })
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Plot potential contours
    contour = ax.contourf(X, Y, potential, 50, cmap='terrain')
    ax.set_title('Potential Energy Surface')
    
    # Add colorbar
    plt.colorbar(contour, ax=ax, label='Potential Energy')
    
    # Plot masses
    ax.add_patch(Circle((0, 0), 0.05, color='blue', label='Mass 1 (k=1.0)'))
    ax.add_patch(Circle((1, 0), 0.02, color='green', label='Mass 2 (k=0.01)'))
    ax.add_patch(Circle((barycenter_x, barycenter_y), 0.01, color='black', label='Barycenter'))
    
    # Plot equilibrium points
    for i, eq in enumerate(equilibrium_points):
        pos = eq['position']
        ax.plot(pos[0], pos[1], marker=eq['marker'], markersize=10, 
                color=eq['color'], label=f"L{i+1}: {eq['stability']}")
    
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.axis('equal')
    ax.grid(True)
    ax.legend()
    
    return fig

# Function to plot a 3D view of the potential surface
def plot_3d_potential():
    system_force, barycenter_x, barycenter_y = setup_system()
    
    # Define the grid
    x = np.linspace(-0.5, 1.5, 100)
    y = np.linspace(-1.0, 1.0, 100)
    X, Y = np.meshgrid(x, y)
    
    # Calculate potential
    potential = np.zeros_like(X)
    
    # Calculate gravitational potential
    r1 = np.sqrt(X**2 + Y**2)
    potential -= 1.0 / r1  # k=1.0
    
    r2 = np.sqrt((X-1)**2 + Y**2)
    potential -= 0.01 / r2  # k=0.01
    
    # Centrifugal potential
    rc = np.sqrt((X-barycenter_x)**2 + (Y-barycenter_y)**2)
    potential += 0.5 * 1.0 * rc**2  # k=1.0
    
    # Find equilibrium points
    initial_guesses = [
        [barycenter_x, 0.0],
        [1.2, 0.0],
        [-0.5, 0.0],
        [0.5, 0.5],
        [0.5, -0.5]
    ]
    
    equilibrium_points = []
    for guess in initial_guesses:
        position, force, converged, _ = newton_raphson(system_force, guess)
        if converged and np.linalg.norm(force) < 1e-6:
            is_duplicate = False
            for p in equilibrium_points:
                if np.linalg.norm(p['position'] - position) < 1e-4:
                    is_duplicate = True
                    break
                    
            if not is_duplicate:
                # Get the potential value at this point
                ix = np.abs(x - position[0]).argmin()
                iy = np.abs(y - position[1]).argmin()
                pot_value = potential[iy, ix]
                
                # Get eigenvalues for stability
                _, jacobian = system_force(position[0], position[1])
                eigenvalues = np.linalg.eigvals(jacobian)
                
                # Classify stability
                if all(np.real(eigenvalues) > 0):
                    stability = "Unstable"
                    marker = 'X'
                    color = 'red'
                elif all(np.real(eigenvalues) < 0):
                    stability = "Stable"
                    marker = 'o'
                    color = 'green'
                else:
                    stability = "Saddle"
                    marker = 's'
                    color = 'orange'
                
                equilibrium_points.append({
                    'position': position,
                    'potential': pot_value,
                    'eigenvalues': eigenvalues,
                    'stability': stability,
                    'marker': marker,
                    'color': color
                })
    
    # Create figure
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    # Plot surface
    # Clip extreme values for better visualization
    vmin = np.percentile(potential, 5)
    vmax = np.percentile(potential, 95)
    surf = ax.plot_surface(X, Y, potential, cmap='terrain', alpha=0.8,
                          vmin=vmin, vmax=vmax)
    
    # Add a color bar
    fig.colorbar(surf, ax=ax, shrink=0.5, aspect=5, label='Potential Energy')
    
    # Plot the equilibrium points
    for i, eq in enumerate(equilibrium_points):
        pos = eq['position']
        ax.scatter(pos[0], pos[1], eq['potential'], 
                  marker=eq['marker'], s=100, color=eq['color'],
                  label=f"L{i+1}: {eq['stability']}")
    
    # Plot masses with vertical lines to the surface
    # Mass 1
    ix1, iy1 = np.abs(x - 0).argmin(), np.abs(y - 0).argmin()
    pot1 = potential[iy1, ix1]
    ax.scatter(0, 0, pot1, s=200, color='blue', label='Mass 1')
    ax.plot([0, 0], [0, 0], [pot1, vmin], 'b--', alpha=0.5)
    
    # Mass 2
    ix2, iy2 = np.abs(x - 1).argmin(), np.abs(y - 0).argmin()
    pot2 = potential[iy2, ix2]
    ax.scatter(1, 0, pot2, s=100, color='green', label='Mass 2')
    ax.plot([1, 1], [0, 0], [pot2, vmin], 'g--', alpha=0.5)
    
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('Potential Energy')
    ax.set_title('3D Potential Energy Surface')
    ax.legend()
    
    return fig

# Main function to generate all plots
def generate_visualizations():
    print("Generating force field visualization...")
    force_fig, points = plot_force_field()
    
    print("\nGenerating potential contour plot...")
    potential_fig = plot_potential()
    
  
    
    
    return force_fig, potential_fig, points

# Run the verification and find equilibrium points
if __name__ == "__main__":
    print("Verifying system against reference values...")
    is_verified = verify_system()
    print(f"System verification: {'Passed' if is_verified else 'Failed'}")
    print("\nFinding equilibrium points...")
    equilibrium_points = find_equilibrium_points()
    print(f"Found {len(equilibrium_points)} equilibrium points")
    generate_visualizations()
