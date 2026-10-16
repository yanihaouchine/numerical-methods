import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import sys
import re
from load_foil import load_foil

def calculate_flow_curves(x, extrados, intrados, num_curves=10):
    """
    Calculate flow curves above and below the airfoil.
    
    Args:
        x: x-coordinates
        extrados: y-coordinates of upper airfoil surface
        intrados: y-coordinates of lower airfoil surface
        num_curves: number of flow curves to generate
        
    Returns:
        Dictionary containing flow curves
    """
    # Calculate min and max heights
    hmin = min(np.min(intrados), np.min(extrados))
    hmax = max(np.max(intrados), np.max(extrados))
    
    # Define boundaries for disturbed airflow
    upper_boundary = 3 * hmax
    lower_boundary = 3 * hmin
    
    # Generate lambda values for flow curves
    lambda_values = np.linspace(0, 1, num_curves)
    
    # Calculate flow curves above the wing (extrados)
    upper_curves = []
    for lam in lambda_values:
        y_curve = (1 - lam) * extrados + lam * upper_boundary
        upper_curves.append(y_curve)
    
    # Calculate flow curves below the wing (intrados)
    lower_curves = []
    for lam in lambda_values:
        y_curve = (1 - lam) * intrados + lam * lower_boundary
        lower_curves.append(y_curve)
    
    return {
        'x': x,
        'upper_curves': upper_curves,
        'lower_curves': lower_curves,
        'lambda_values': lambda_values,
        'upper_boundary': upper_boundary,
        'lower_boundary': lower_boundary
    }

def calculate_curve_lengths(x, curves):
    """
    Calculate the length of each flow curve.
    
    Args:
        x: x-coordinates
        curves: list of y-coordinates for curves
        
    Returns:
        Array of curve lengths
    """
    lengths = []
    
    for curve in curves:
        # Make sure we have the same length for x and curve
        if len(x) != len(curve):
            raise ValueError(f"Mismatch in array lengths: x={len(x)}, curve={len(curve)}")
            
        # Calculate length using segment-by-segment approach
        total_length = 0
        for i in range(len(x) - 1):  # Stop at len(x)-1 to avoid index error
            dx = x[i+1] - x[i]
            dy = curve[i+1] - curve[i]
            segment_length = np.sqrt(dx**2 + dy**2)
            total_length += segment_length
        
        lengths.append(total_length)
    
    return np.array(lengths)

def calculate_velocities(curve_lengths, reference_velocity=1.0):
    curve_lengths = np.asarray(curve_lengths, dtype=float)

    reference_length = curve_lengths.max()
    velocities = reference_velocity * (reference_length / curve_lengths)

    return velocities          # ndarray 1 D, plus jamais de liste


def calculate_pressure_coefficients(velocities, reference_velocity=1.0):
    velocities = np.asarray(velocities, dtype=float)
    return [1 - (velocities / reference_velocity)**2]
  





def calculate_lift(x, upper_pressure_coeffs, lower_pressure_coeffs):
    """
    Calcule le coefficient de portance C_L.

    x                  : abscisses le long de l'extrados (N points)
    upper_pressure_coeffs[0] : Cp sur l'extrados   (M points ou scalaire)
    lower_pressure_coeffs[0] : Cp sur l'intrados   (M points ou scalaire)
    Retour : C_L
    """
    x = np.asarray(x, dtype=float)

    # --- 1. On récupère le Cp "surface" (λ = 0) ------------------------
    Cp_up  = np.asarray(upper_pressure_coeffs[0], dtype=float).ravel()
    Cp_low = np.asarray(lower_pressure_coeffs[0], dtype=float).ravel()

    # --- 2. On ajuste chaque Cp à la longueur de x ---------------------
    def stretch(Cp):
        if Cp.size == 1:                      # scalaire → on duplique
            return np.full_like(x, Cp.item())
        if Cp.size == x.size:                # déjà à la bonne taille
            return Cp
        # sinon : interpolation linéaire
        s = np.linspace(x.min(), x.max(), Cp.size)
        return np.interp(x, s, Cp)

    Cp_up  = stretch(Cp_up)
    Cp_low = stretch(Cp_low)

    # --- 3. Intégration trapézoïdale -----------------------------------
    dCp = Cp_low - Cp_up
    cl  = np.trapz(dCp, x) / (x[-1] - x[0])   # normalisation par la corde
    return cl


def visualize_airflow(x_extrados, y_extrados,
                      x_intrados, y_intrados,
                      flow_data, velocities):

    import matplotlib.pyplot as plt
    import numpy as np

    # ------------------------------------------------------------------ #
    # 1. Figure et Axe explicites                                        #
    # ------------------------------------------------------------------ #
    fig, ax = plt.subplots(figsize=(14, 8))

    # --- Interpolation éventuelle de l’intrados ----------------------- #
    if y_intrados.size != x_extrados.size:
        y_intrados = np.interp(x_extrados, x_intrados, y_intrados)

    # --- Profil plein gris ------------------------------------------- #
    x_poly = np.concatenate((x_extrados, x_extrados[::-1]))
    y_poly = np.concatenate((y_extrados, y_intrados[::-1]))
    ax.fill(x_poly, y_poly, color='lightgray', alpha=0.8, label='Airfoil')
    ax.plot(x_extrados, y_extrados, 'k-', lw=2)
    ax.plot(x_extrados, y_intrados, 'k-', lw=2)

    # ------------------------------------------------------------------ #
    # 2. Lignes de courant colorées                                     #
    # ------------------------------------------------------------------ #
    cmap = plt.cm.jet
    up_vel  = velocities['upper_velocities']
    low_vel = velocities['lower_velocities']
    n_up  = plt.Normalize(min(up_vel),  max(up_vel))
    n_low = plt.Normalize(min(low_vel), max(low_vel))

    for i, curve in enumerate(flow_data['upper_curves']):
        ax.plot(x_extrados, curve, color=cmap(n_up(up_vel[i])), alpha=0.7)
    for i, curve in enumerate(flow_data['lower_curves']):
        ax.plot(x_intrados, curve, color=cmap(n_low(low_vel[i])), alpha=0.7)

    # ------------------------------------------------------------------ #
    # 3. Barre de couleurs explicite                                     #
    # ------------------------------------------------------------------ #
    sm = plt.cm.ScalarMappable(
            cmap=cmap,
            norm=plt.Normalize(min(min(up_vel), min(low_vel)),
                               max(max(up_vel), max(low_vel))))
    sm.set_array([])
    fig.colorbar(sm, ax=ax, label='Relative Air Velocity')   # <- ICI

    # ------------------------------------------------------------------ #
    # 4. Finitions                                                       #
    # ------------------------------------------------------------------ #
    ax.set_title('Airflow Around the Airfoil')
    ax.set_xlabel('x/c')
    ax.set_ylabel('y/c')
    ax.set_aspect('equal')
    ax.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.show()



def visualize_pressure(x_extrados, y_extrados,
                       x_intrados, y_intrados,
                       upper_pressure, lower_pressure):
    """
    Affiche le profil + la distribution de Cp.
    Les longueurs des tableaux extrados/intrados peuvent différer :
    on interpole l’intrados sur la même grille x que l’extrados.
    """
    import numpy as np
    import matplotlib.pyplot as plt

    # ------------------------------------------------------------------ #
    # 1. Mise à la même grille                                           #
    # ------------------------------------------------------------------ #
    if y_intrados.size != x_extrados.size:
        y_intrados = np.interp(x_extrados, x_intrados, y_intrados)

    # idem pour les coefficients de pression (si ≠ 33)
    def stretch(Cp):
        Cp = np.asarray(Cp[0], dtype=float).ravel()  # on ne garde que λ=0
        if Cp.size == 1:
            return np.full_like(x_extrados, Cp.item())
        if Cp.size == x_extrados.size:
            return Cp
        s = np.linspace(x_extrados.min(), x_extrados.max(), Cp.size)
        return np.interp(x_extrados, s, Cp)

    Cp_up  = stretch(upper_pressure)
    Cp_low = stretch(lower_pressure)

    # ------------------------------------------------------------------ #
    # 2. Figure avec deux sous‑graphiques                                #
    # ------------------------------------------------------------------ #
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8),
                                   gridspec_kw={'height_ratios': [1, 1.3]})

    # --- Profil plein gris ------------------------------------------- #
    x_poly = np.concatenate((x_extrados, x_extrados[::-1]))
    y_poly = np.concatenate((y_extrados, y_intrados[::-1]))
    ax1.fill(x_poly, y_poly, color='lightgray', alpha=0.8)
    ax1.plot(x_extrados, y_extrados, 'k-', lw=1.5)
    ax1.plot(x_extrados, y_intrados, 'k-', lw=1.5)
    ax1.set_aspect('equal')
    ax1.set_title('Airfoil Profile')
    ax1.set_xlabel('x/c')
    ax1.set_ylabel('y/c')
    ax1.grid(True, ls='--', alpha=0.7)

    # --- Distribution de Cp ------------------------------------------ #
    ax2.plot(x_extrados, Cp_up,  'r-', label='Extrados')
    ax2.plot(x_extrados, Cp_low, 'b-', label='Intrados')
    ax2.fill_between(x_extrados, Cp_up, Cp_low,
                     color='green', alpha=0.3,
                     label='ΔCp (portance)')
    ax2.invert_yaxis()  # convention aéronautique
    ax2.set_title('Pressure Coefficient Distribution')
    ax2.set_xlabel('x/c')
    ax2.set_ylabel('Cp')
    ax2.grid(True, ls='--', alpha=0.7)
    ax2.legend()

    plt.tight_layout()
    plt.show()


def create_pressure_contour_plot(x, extrados, intrados, flow_data, upper_pressure, lower_pressure):
    """
    Create a contour plot of pressure distribution around the airfoil.
    
    Args:
        x: x-coordinates for extrados
        extrados: upper airfoil surface
        intrados: lower airfoil surface (might have different length than extrados)
        flow_data: dictionary with flow curve data
        upper_pressure: pressure coefficients on upper curves
        lower_pressure: pressure coefficients on lower curves
    """
    plt.figure(figsize=(14, 8))
    
    # First, make sure intrados is interpolated to the same x-grid as extrados
    if len(intrados) != len(x):
        # We need the x-coordinates for intrados
        ix = flow_data.get('x', np.linspace(0, 1, len(intrados)))
        # Interpolate intrados to match x-grid of extrados
        intrados_interp = np.interp(x, ix, intrados)
    else:
        intrados_interp = intrados
    
    # Generate a grid for pressure visualization (fewer points for efficiency)
    x_grid = np.linspace(min(x), max(x), 50)
    
    # Create empty grid for pressure data
    y_min = flow_data['lower_boundary'] - 0.05
    y_max = flow_data['upper_boundary'] + 0.05
    y_grid = np.linspace(y_min, y_max, 50)
    pressure_grid = np.zeros((len(y_grid), len(x_grid)))
    
    # Convert pressure coefficients to numpy arrays if they aren't already
    upper_pressure_array = np.asarray(upper_pressure[0]).flatten()
    lower_pressure_array = np.asarray(lower_pressure[0]).flatten()
    
    # Ensure we have pressure coefficients for each x-coordinate
    if len(upper_pressure_array) == 1:
        upper_pressure_array = np.full(len(x), upper_pressure_array[0])
    if len(lower_pressure_array) == 1:
        lower_pressure_array = np.full(len(x), lower_pressure_array[0])
    
    # Interpolate pressure coefficients to match x-grid if needed
    if len(upper_pressure_array) != len(x):
        x_orig = np.linspace(min(x), max(x), len(upper_pressure_array))
        upper_pressure_array = np.interp(x, x_orig, upper_pressure_array)
    if len(lower_pressure_array) != len(x):
        x_orig = np.linspace(min(x), max(x), len(lower_pressure_array))
        lower_pressure_array = np.interp(x, x_orig, lower_pressure_array)
    
    # Interpolate extrados and intrados to match x_grid
    extrados_grid = np.interp(x_grid, x, extrados)
    intrados_grid = np.interp(x_grid, x, intrados_interp)
    
    # Create simplified airfoil mask
    airfoil_mask = np.zeros((len(y_grid), len(x_grid)), dtype=bool)
    
    # Fill pressure grid with a simpler approach
    X_grid, Y_grid = np.meshgrid(x_grid, y_grid)
    
    # For each point in the grid
    for i in range(len(y_grid)):
        for j in range(len(x_grid)):
            x_val = x_grid[j]
            y_val = y_grid[i]
            
            # Determine if this point is inside the airfoil
            if intrados_grid[j] < y_val < extrados_grid[j]:
                airfoil_mask[i, j] = True
                pressure_grid[i, j] = np.nan
                continue
            
            # Simple distance-based interpolation for pressure
            # Above the airfoil
            if y_val >= extrados_grid[j]:
                dist_from_surface = y_val - extrados_grid[j]
                total_dist = flow_data['upper_boundary'] - extrados_grid[j]
                factor = min(dist_from_surface / total_dist, 1.0) if total_dist > 0 else 1.0
                # Blend pressure coefficient (higher pressure far from airfoil)
                pressure_grid[i, j] = upper_pressure_array[0] * (1 - factor) + 0 * factor
            
            # Below the airfoil
            elif y_val <= intrados_grid[j]:
                dist_from_surface = intrados_grid[j] - y_val
                total_dist = intrados_grid[j] - flow_data['lower_boundary']
                factor = min(dist_from_surface / total_dist, 1.0) if total_dist > 0 else 1.0
                # Blend pressure coefficient (higher pressure far from airfoil)
                pressure_grid[i, j] = lower_pressure_array[0] * (1 - factor) + 0 * factor
    
    # Apply mask to pressure grid
    masked_pressure = np.ma.array(pressure_grid, mask=airfoil_mask)
    
    # Plot contour
    contour = plt.contourf(X_grid, Y_grid, masked_pressure, 20, cmap='jet')
    plt.colorbar(contour, label='Pressure Coefficient (Cp)')
    
    # Plot airfoil
    plt.fill(x, extrados, x[::-1], intrados_interp[::-1], color='gray', alpha=1.0)
    plt.plot(x, extrados, 'k-', linewidth=1.5)
    plt.plot(x, intrados_interp, 'k-', linewidth=1.5)
    
    plt.title('Pressure Distribution Around Airfoil')
    plt.xlabel('x/c (chord position)')
    plt.ylabel('y/c (relative position)')
    plt.axis('equal')
    plt.grid(True, linestyle='--', alpha=0.3)
    plt.xlim(-0.05, 1.05)
    plt.ylim(y_min, y_max)
    plt.tight_layout()
    plt.show()

def analyze_angle_of_attack(filename, angles=[-5, 0, 5, 10, 15]):
    """
    Analyze the effect of angle of attack on lift coefficient.
    
    Args:
        filename: path to airfoil data file
        angles: list of angles of attack to analyze (in degrees)
    """
    # Load airfoil data
    (dim, ex, ey, ix, iy) = load_foil(filename)
    
    # Calculate chord line and chord length
    chord_start = [ex[0], ey[0]]  # Leading edge
    chord_end = [ex[-1], ey[-1]]  # Trailing edge
    chord_length = np.sqrt((chord_end[0] - chord_start[0])**2 + (chord_end[1] - chord_start[1])**2)
    
    # Initialize lists to store results
    lift_coefficients = []
    
    for angle in angles:
        # Convert angle to radians
        angle_rad = np.radians(angle)
        
        # Rotate airfoil coordinates
        rotated_ex = []
        rotated_ey = []
        rotated_ix = []
        rotated_iy = []
        
        for i in range(len(ex)):
            # Translate to origin (leading edge)
            x_translated = ex[i] - chord_start[0]
            y_translated = ey[i] - chord_start[1]
            
            # Rotate
            x_rotated = x_translated * np.cos(angle_rad) - y_translated * np.sin(angle_rad)
            y_rotated = x_translated * np.sin(angle_rad) + y_translated * np.cos(angle_rad)
            
            # Translate back
            rotated_ex.append(x_rotated + chord_start[0])
            rotated_ey.append(y_rotated + chord_start[1])
        
        for i in range(len(ix)):
            # Translate to origin (leading edge)
            x_translated = ix[i] - chord_start[0]
            y_translated = iy[i] - chord_start[1]
            
            # Rotate
            x_rotated = x_translated * np.cos(angle_rad) - y_translated * np.sin(angle_rad)
            y_rotated = x_translated * np.sin(angle_rad) + y_translated * np.cos(angle_rad)
            
            # Translate back
            rotated_ix.append(x_rotated + chord_start[0])
            rotated_iy.append(y_rotated + chord_start[1])
        
        # Convert to numpy array
        rotated_ex = np.array(rotated_ex)
        rotated_ey = np.array(rotated_ey)
        rotated_ix = np.array(rotated_ix)
        rotated_iy = np.array(rotated_iy)
        
        # Calculate flow curves
        flow_data = calculate_flow_curves(rotated_ex, rotated_ey, rotated_iy)
        
        # Calculate curve lengths
        upper_lengths = calculate_curve_lengths(rotated_ex, flow_data['upper_curves'])
        lower_lengths = calculate_curve_lengths(rotated_ix, flow_data['lower_curves'])
        
        # Calculate velocities
        upper_velocities = calculate_velocities(upper_lengths)
        lower_velocities = calculate_velocities(lower_lengths)
        
        # Calculate pressure coefficients
        upper_pressure_coeffs = calculate_pressure_coefficients(upper_velocities)
        lower_pressure_coeffs = calculate_pressure_coefficients(lower_velocities)
        
        # Calculate lift coefficient
        lift_coefficient = calculate_lift(rotated_ex, upper_pressure_coeffs, lower_pressure_coeffs)
        lift_coefficients.append(lift_coefficient)
        
        print(f"Angle of Attack: {angle}° - Lift Coefficient: {lift_coefficient:.4f}")
    
    # Plot lift coefficient vs angle of attack
    plt.figure(figsize=(10, 6))
    plt.plot(angles, lift_coefficients, 'bo-', linewidth=2)
    plt.xlabel('Angle of Attack (degrees)')
    plt.ylabel('Lift Coefficient')
    plt.title('Lift Coefficient vs Angle of Attack')
    plt.grid(True)
    plt.tight_layout()
    plt.show()
    
    return {'angles': angles, 'lift_coefficients': lift_coefficients}


    
def analyze_airfoil_pressure(filename):
    """
    Complete analysis of airfoil pressure distribution.
    
    Args:
        filename: path to airfoil data file
    """
    # Load airfoil data
    print(f"Loading airfoil data from {filename}...")
    (dim, ex, ey, ix, iy) = load_foil(filename)
    
    print(f"Airfoil data loaded. Extrados points: {len(ex)}, Intrados points: {len(ix)}")
    
    # Calculate flow curves
    print("Calculating flow curves...")
    flow_data = calculate_flow_curves(ex, ey, iy)
    
    # Debug info
    print(f"Upper curves: {len(flow_data['upper_curves'])}, points per curve: {len(flow_data['upper_curves'][0])}")
    print(f"Lower curves: {len(flow_data['lower_curves'])}, points per curve: {len(flow_data['lower_curves'][0])}")
    print(f"x points: {len(flow_data['x'])}")
    
    # Calculate curve lengths
    print("Calculating curve lengths...")
    upper_lengths = calculate_curve_lengths(ex, flow_data['upper_curves'])
    lower_lengths = calculate_curve_lengths(ix, flow_data['lower_curves'])
    
    print(f"Upper lengths shape: {np.shape(upper_lengths)}")
    print(f"Lower lengths shape: {np.shape(lower_lengths)}")
    
    # Calculate velocities
    print("Calculating velocities...")
    upper_velocities = calculate_velocities(upper_lengths)
    lower_velocities = calculate_velocities(lower_lengths)
    
    print(f"Upper velocities type: {type(upper_velocities)}")
    print(f"Upper velocities shape: {np.shape(upper_velocities) if hasattr(upper_velocities, 'shape') else len(upper_velocities)}")
    
    velocities = {
        'upper_velocities': upper_velocities,
        'lower_velocities': lower_velocities
    }
    
    # Calculate pressure coefficients
    print("Calculating pressure coefficients...")
    upper_pressure_coeffs = calculate_pressure_coefficients(upper_velocities)
    lower_pressure_coeffs = calculate_pressure_coefficients(lower_velocities)
    
    print(f"Upper pressure coeffs type: {type(upper_pressure_coeffs)}")
    print(f"Upper pressure coeffs shape: {np.shape(upper_pressure_coeffs) if hasattr(upper_pressure_coeffs, 'shape') else len(upper_pressure_coeffs)}")
    
    # Calculate lift coefficient
    print("Calculating lift coefficient...")
    lift_coefficient = calculate_lift(ex, upper_pressure_coeffs, lower_pressure_coeffs)
    
    # Print results
    print(f"\nAirfoil Analysis Results:")
    print(f"------------------------")
    print(f"Min/Max Height: {min(np.min(iy), np.min(ey)):.4f} / {max(np.max(iy), np.max(ey)):.4f}")
    print(f"Lift Coefficient: {lift_coefficient:.4f}")
    
    if hasattr(upper_velocities, "__len__") and len(upper_velocities) > 0:
        if hasattr(upper_velocities[0], "__len__"):
            print(f"Max Velocity (Upper Surface): {max(upper_velocities[0]):.4f}")
        else:
            print(f"Max Velocity (Upper Surface): {max(upper_velocities):.4f}")
    
    if hasattr(lower_velocities, "__len__") and len(lower_velocities) > 0:
        if hasattr(lower_velocities[0], "__len__"):
            print(f"Max Velocity (Lower Surface): {max(lower_velocities[0]):.4f}")
        else:
            print(f"Max Velocity (Lower Surface): {max(lower_velocities):.4f}")
    
    # Visualize airflow
    print("\nVisualizing airflow...")
    visualize_airflow(ex, ey, ix, iy, flow_data, velocities)
    
    # Visualize pressure distribution
    print("Visualizing pressure distribution...")
    visualize_pressure(ex, ey, ix, iy, upper_pressure_coeffs, lower_pressure_coeffs)
    
    return {
        'lift': lift_coefficient,
        'upper_pressure': upper_pressure_coeffs,
        'lower_pressure': lower_pressure_coeffs,
        'upper_velocity': upper_velocities,
        'lower_velocity': lower_velocities,
        'flow_data': flow_data
    }
if __name__ == "__main__":
    if len(sys.argv) == 1:
        print("Please provide an input file")
        exit(1)
    
    try:
        # Complete airfoil analysis
        results = analyze_airfoil_pressure(sys.argv[1])
        
        # Get airfoil data for additional visualizations
        (dim, ex, ey, ix, iy) = load_foil(sys.argv[1])
        
        # # Create pressure contour plot
        # print("\nCreating pressure contour plot...")
        # create_pressure_contour_plot(ex, ey, iy, results['flow_data'], results['upper_pressure'], results['lower_pressure'])
        
        # # Analyze effect of angle of attack
        # print("\nAnalyzing effect of angle of attack...")
        # analyze_angle_of_attack(sys.argv[1])
        
    except Exception as e:
        print(f"Error during analysis: {e}")
        import traceback
        traceback.print_exc()
