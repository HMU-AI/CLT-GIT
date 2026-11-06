

import os
import numpy as np
import pandas as pd
import cv2
from tqdm import tqdm


INPUT_FOLDER = '../expr_csv_aligned'

OUTPUT_ROOT_DIR = '../output_heatmaps'

os.makedirs(OUTPUT_ROOT_DIR, exist_ok=True)


def genes_to_image(expression_vector, genes_per_row=136):
    """
    expression_vector: (N,) numpy array
    return: uint8 heatmap (H, W)
    """
    n_genes = len(expression_vector)
    n_rows = int(np.ceil(n_genes / genes_per_row))

    # padding
    padded = np.pad(
        expression_vector,
        (0, n_rows * genes_per_row - n_genes),
        mode='constant'
    )

    image_2d = padded.reshape(n_rows, genes_per_row)


    low, high = np.quantile(image_2d, [0.1, 0.9])
    image_2d = np.clip(image_2d, low, high)


    image_2d = (image_2d - low) / (high - low + 1e-8)
    image_2d = (image_2d * 255).astype(np.uint8)

    return image_2d


for filename in os.listdir(INPUT_FOLDER):

    if not filename.endswith('.csv'):
        continue
    

    csv_path = os.path.join(INPUT_FOLDER, filename)

    csv_name = os.path.splitext(filename)[0]
    output_dir = os.path.join(OUTPUT_ROOT_DIR, csv_name)
    os.makedirs(output_dir, exist_ok=True)

    print(f"\nprocessing: {filename}")
    print(f"output root: {output_dir}")


    df = pd.read_csv(csv_path, index_col=0, low_memory=False)
    df = df.T
    df.columns = df.columns.astype(str).str.upper()
    norm_df = df.copy()

    print(f"size: {norm_df.shape}")


    for sample_name in tqdm(norm_df.index, desc=f" {csv_name}"):
        vec = norm_df.loc[sample_name].values


        img_gray = genes_to_image(vec, genes_per_row=136)


        img_color = cv2.applyColorMap(img_gray, cv2.COLORMAP_VIRIDIS)


        safe_name = str(sample_name).replace('/', '_').replace(':', '_')
        save_path = os.path.join(output_dir, f"{safe_name}_heatmap.png")


        cv2.imwrite(save_path, img_color)

    print(f"✔  {csv_name}： {len(norm_df)} ")

print(f"\n🎉 ok:{OUTPUT_ROOT_DIR}")