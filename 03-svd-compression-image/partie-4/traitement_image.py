import numpy as np
from PIL import Image
#from sklearn.utils.extmath import randomized_svd
import matplotlib.pyplot as plt
import os

def randomized_svd(A, n_components, n_oversamples=10, n_iter=4, random_state=None):
    """
    Computes a randomized SVD on the input matrix A.
    
    Parameters:
    -----------
    A : array-like, shape (m, n)
        Input matrix to decompose
    n_components : int
        Number of singular values and vectors to extract
    n_oversamples : int, optional (default=10)
        Additional number of random vectors to sample the range of A
    n_iter : int, optional (default=4)
        Number of power iterations
    random_state : int, RandomState instance or None, optional (default=None)
        The seed of the pseudo random number generator
        
    Returns:
    --------
    U : array-like, shape (m, n_components)
        Left singular vectors
    s : array-like, shape (n_components,)
        Singular values
    Vt : array-like, shape (n_components, n)
        Right singular vectors (transposed)
    """
    random_state = np.random.RandomState(random_state)
    
    m, n = A.shape
    
    if n_components > min(m, n):
        raise ValueError("n_components must be <= min(m, n)")
    
    # Step 1: Sample the range of A using a random matrix
    n_random = n_components + n_oversamples
    Q = random_state.normal(size=(n, n_random))
    
    # Step 2: Project A onto the random vectors Q
    Y = A @ Q
    
    # Step 3: Apply power iterations to increase accuracy
    for _ in range(n_iter):
        Y = A @ (A.T @ Y)
    
    # Step 4: Orthogonalize Y using QR decomposition
    Q, _ = linalg.qr(Y, mode='economic')
    
    # Step 5: Project matrix A onto the orthogonal basis Q
    B = Q.T @ A
    
    # Step 6: Compute the SVD of the smaller matrix B
    Uhat, s, Vt = linalg.svd(B, full_matrices=False)
    
    # Step 7: Transform Uhat to get the left singular vectors of A
    U = Q @ Uhat
    
    # Return only the requested number of components
    return U[:, :n_components], s[:n_components], Vt[:n_components, :]


 
    
def compress_image_partial_svd_rgb(image_path, k, output_file):
    """
    Charge une image RGB, applique la SVD partielle sur chaque canal (R, G, B)
    et sauvegarde la décomposition (matrices U, s, Vt) ainsi que les dimensions
    de l'image dans un fichier .npz.
    
    Paramètres :
      - image_path  : chemin vers l'image d'entrée (format RGB)
      - k           : rang de la décomposition SVD
      - output_file : chemin du fichier de sortie (.npz)
    """
    # Lecture de l'image et conversion en tableau de floats
    img = Image.open(image_path).convert('RGB')
    arr = np.array(img, dtype=np.float64)  # shape (H, W, 3)
    H, W, _ = arr.shape
    
    # Séparation des canaux R, G, B
    R_channel = arr[..., 0]
    G_channel = arr[..., 1]
    B_channel = arr[..., 2]
    
    # Fonction pour effectuer la SVD partielle sur un canal
    def partial_svd_channel(channel, k):
        U, s, Vt = randomized_svd(channel, n_components=k)
        return U, s, Vt
    
    # Application de la SVD partielle sur chaque canal
    R_U, R_s, R_Vt = partial_svd_channel(R_channel, k)
    G_U, G_s, G_Vt = partial_svd_channel(G_channel, k)
    B_U, B_s, B_Vt = partial_svd_channel(B_channel, k)
    
    # Sauvegarde des décompositions et de la taille de l'image dans un fichier .npz compressé
    np.savez_compressed(
        output_file,
        H=H, W=W,
        R_U=R_U, R_s=R_s, R_Vt=R_Vt,
        G_U=G_U, G_s=G_s, G_Vt=G_Vt,
        B_U=B_U, B_s=B_s, B_Vt=B_Vt
    )
    print(f"[COMPRESSION] k={k} -> Données sauvegardées dans {output_file}.")

def decompress_image_partial_svd_rgb(input_file, output_image):
    """
    Recharge la décomposition SVD partielle stockée dans input_file (.npz),
    reconstruit chaque canal (R, G, B), recompose l'image RGB et sauvegarde le résultat.
    
    Paramètres :
      - input_file  : chemin vers le fichier .npz contenant la décomposition
      - output_image: chemin du fichier de sortie (image reconstruite)
    """
    data = np.load(input_file)
    H = int(data["H"])
    W = int(data["W"])
    
    # Fonction de reconstruction d'un canal à partir de (U, s, Vt)
    def reconstruct_channel(U, s, Vt):
        return U @ np.diag(s) @ Vt
    
    # Reconstruction de chaque canal
    R = reconstruct_channel(data["R_U"], data["R_s"], data["R_Vt"])
    G = reconstruct_channel(data["G_U"], data["G_s"], data["G_Vt"])
    B = reconstruct_channel(data["B_U"], data["B_s"], data["B_Vt"])
    
    # Recomposition de l'image RGB
    reconstructed = np.zeros((H, W, 3), dtype=np.float64)
    reconstructed[..., 0] = R
    reconstructed[..., 1] = G
    reconstructed[..., 2] = B
    
    # Recadrage des valeurs dans [0, 255] et conversion en uint8 pour l'affichage
    reconstructed = np.clip(reconstructed, 0, 255).astype(np.uint8)
    
    # Sauvegarde de l'image reconstruite
    img_out = Image.fromarray(reconstructed, mode='RGB')
    img_out.save(output_image)
    print(f"[DECOMPRESSION] Image sauvegardée dans {output_image}.")

def mse(image1, image2):
    """
    Calcul de l'erreur quadratique moyenne (MSE) entre deux images.
    Les images doivent être de mêmes dimensions et de type float ou int.
    """
    return np.mean((image1 - image2) ** 2)

def psnr(original, reconstructed):
    """
    Calcul du Peak Signal-to-Noise Ratio (PSNR).
    """
    mse_val = mse(original, reconstructed)
    if mse_val == 0:
        return float('inf')
    max_val = 255.0  # si vos images sont dans [0, 255]
    return 20 * np.log10(max_val / np.sqrt(mse_val))


if __name__ == "__main__":
    # -------------------------------------------------------------------------
    # Exemple d'utilisation simple : compresser, puis décompresser une image
    # -------------------------------------------------------------------------
    image_path = os.path.join(os.path.dirname(__file__), "..", "img", "p3_takeoff_base.png")  # Chemin vers l'image d'entrée (RGB)
    output_npz = "compressed_svd.npz"
    output_image = "decompressed_svd.png"
    k_example = 5
    
    compress_image_partial_svd_rgb(image_path, k_example, output_npz)
    decompress_image_partial_svd_rgb(output_npz, output_image)
    img_orig = Image.open(image_path).convert('RGB')
    img_recon = Image.open(output_image)
    
    plt.figure(figsize=(12, 6))
    plt.subplot(1, 2, 1)
    plt.imshow(img_orig)
    plt.title("Image originale")
    plt.axis('off')
    
    plt.subplot(1, 2, 2)
    plt.imshow(img_recon)
    plt.title(f"Image reconstruite (k_example={k_example})")
    plt.axis('off')
    
    plt.tight_layout()
    plt.show()

    # -------------------------------------------------------------------------
    # Évaluation de la qualité de compression (MSE) pour différentes valeurs de k
    # -------------------------------------------------------------------------
    original_img = Image.open(image_path).convert('RGB')
    original_arr = np.array(original_img, dtype=np.float64)

    ks = range(1, 51, 1)  
    errors = []
    
    # Dossier temporaire pour stocker les fichiers intermédiaires
    temp_dir = "temp_svd"
    os.makedirs(temp_dir, exist_ok=True)

    for k in ks:
        # Fichiers .npz et .png temporaires
        tmp_npz = os.path.join(temp_dir, f"temp_{k}.npz")
        tmp_png = os.path.join(temp_dir, f"temp_{k}.png")
        
        # Compression
        compress_image_partial_svd_rgb(image_path, k, tmp_npz)
        # Décompression
        decompress_image_partial_svd_rgb(tmp_npz, tmp_png)
        
        # Lecture de l'image reconstruite
        reconstructed_img = Image.open(tmp_png).convert('RGB')
        reconstructed_arr = np.array(reconstructed_img, dtype=np.float64)
        
        # Calcul de la MSE
        error = mse(original_arr, reconstructed_arr)
        errors.append(error)
    
    # Tracé de la MSE en fonction de k
    plt.figure(figsize=(8, 5))
    plt.plot(ks, errors, marker='o')
    plt.xlabel("Rang k")
    plt.ylabel("Erreur MSE")
    plt.title("Effet du rang k sur la qualité de la compression (MSE)")
    plt.grid(True)
    plt.show()

    # (Optionnel) Affichage final de la PSNR pour le dernier k
    final_psnr = psnr(original_arr, reconstructed_arr)
    print(f"PSNR pour k={k}: {final_psnr:.2f} dB")

    # -------------------------------------------------------------------------
    # Nettoyage éventuel des fichiers temporaires (optionnel)
    # -------------------------------------------------------------------------
    # import shutil
    # shutil.rmtree(temp_dir)  # Supprime le dossier et tout son contenu
