import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from umap import UMAP
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import numpy as np

def run_pca(adata, n_components=50, standardize=True, copy=False, key_added='X_pca'):
    """
    Compute PCA and store results in `adata.obsm[key_added]` and `adata.uns['pca']`.

    Parameters:
        adata (AnnData): AnnData object.
        n_components (int): Number of PCA components to compute.
        standardize (bool): Whether to standardize features before PCA.
        copy (bool): Whether to operate on a copy.
        key_added (str): Key to store PCA coordinates in `adata.obsm`.

    Returns:
        AnnData: Modified AnnData object (if copy=True, a copy).
    """
    if copy:
        adata = adata.copy()

    X = adata.X.copy()
    X = np.nan_to_num(X, nan=0, posinf=0, neginf=0)

    if standardize:
        scaler = StandardScaler()
        X = scaler.fit_transform(X)

    pca = PCA(n_components=n_components)
    X_pca = pca.fit_transform(X)

    adata.obsm[key_added] = X_pca
    adata.uns['pca'] = {
        'variance_ratio': pca.explained_variance_ratio_,
        'variance': pca.explained_variance_,
        'total_variance': np.sum(pca.explained_variance_ratio_),
        'components': pca.components_,
        'means': pca.mean_,
    }

    return adata

def plot_pca(
    adata,
    obs_column: str,
    obs_values: list = None,
    var_names: list = None,
    n_components: int = 2,
    standardize: bool = True,
    title: str = None,
    show: bool = True
):
    """
    Filter AnnData and plot PCA colored by obs_column.

    Parameters:
        adata (AnnData): The input data.
        obs_column (str): Column in adata.obs to color by.
        obs_values (list): Values to include from obs_column. If None, use all.
        var_names (list): Genes/features to include. If None, use all.
        n_components (int): Number of PCA components to compute.
        standardize (bool): Standardize features before PCA.
        title (str): Optional plot title.
        show (bool): Whether to display the plot.

    Returns:
        None
    """
    # Subset by obs
    if obs_values is not None:
        adata = adata[adata.obs[obs_column].isin(obs_values)].copy()

    # Subset by features
    if var_names is not None:
        adata = adata[:, var_names].copy()

    # Run PCA
    adata = run_pca(adata, n_components=n_components, standardize=standardize)

    # Prepare dataframe
    df = pd.DataFrame(adata.obsm["X_pca"][:, :2], columns=["PC1", "PC2"])
    df[obs_column] = adata.obs[obs_column].values

    # Plot
    plt.figure(figsize=(8, 6))
    sns.scatterplot(
        data=df,
        x="PC1",
        y="PC2",
        hue=obs_column,
        palette="Set2",
        edgecolor="black",
        alpha=0.75
    )

    var_ratio = adata.uns['pca']['variance_ratio']
    plt.xlabel(f"PC1 ({var_ratio[0]*100:.1f}%)")
    plt.ylabel(f"PC2 ({var_ratio[1]*100:.1f}%)")
    plt.title(title or f"PCA colored by {obs_column}")
    plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left", title=obs_column)
    plt.tight_layout()

    if show:
        plt.show()

def run_umap(adata, n_components=2, standardize=True, copy=False, key_added='X_umap'):
    """
    Compute UMAP and store in `adata.obsm[key_added]`.

    Parameters:
        adata (AnnData): AnnData object.
        n_components (int): Dimensionality of the UMAP projection.
        standardize (bool): Standardize features before UMAP.
        copy (bool): Return a new AnnData object instead of modifying in-place.
        key_added (str): Key for `adata.obsm` to store the result.

    Returns:
        AnnData: Modified or new AnnData object with UMAP in obsm.
    """
    if copy:
        adata = adata.copy()

    X = adata.X.copy()
    X = np.nan_to_num(X, nan=0, posinf=0, neginf=0)

    if standardize:
        scaler = StandardScaler()
        X = scaler.fit_transform(X)

    reducer = UMAP(n_components=n_components, random_state=42)
    X_umap = reducer.fit_transform(X)

    adata.obsm[key_added] = X_umap
    adata.uns['umap'] = {
        'n_components': n_components,
        'random_state': 42
    }

    return adata

def plot_umap(
    adata,
    obs_column: str,
    obs_values: list = None,
    var_names: list = None,
    n_components: int = 2,
    standardize: bool = True,
    title: str = None,
    show: bool = True
):
    """
    Filter AnnData and plot UMAP colored by obs_column.

    Parameters:
        adata (AnnData): Input AnnData object.
        obs_column (str): Column in adata.obs to color by.
        obs_values (list): List of values to include from obs_column. If None, include all.
        var_names (list): List of gene/features to include. If None, use all.
        n_components (int): Dimensionality of UMAP.
        standardize (bool): Whether to standardize the data before UMAP.
        title (str): Optional title for plot.
        show (bool): Whether to display the plot.

    Returns:
        None
    """
    # Subset by obs
    if obs_values is not None:
        adata = adata[adata.obs[obs_column].isin(obs_values)].copy()

    # Subset by var
    if var_names is not None:
        adata = adata[:, var_names].copy()

    # Run UMAP
    adata = run_umap(adata, n_components=n_components, standardize=standardize)

    # Prepare dataframe
    coords = adata.obsm["X_umap"][:, :2]
    df = pd.DataFrame(coords, columns=["UMAP1", "UMAP2"])
    df[obs_column] = adata.obs[obs_column].values

    # Plot
    plt.figure(figsize=(8, 6))
    sns.scatterplot(
        data=df,
        x="UMAP1",
        y="UMAP2",
        hue=obs_column,
        palette="Set2",
        edgecolor="black",
        alpha=0.75
    )
    plt.title(title or f"UMAP colored by {obs_column}")
    plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left", title=obs_column)
    plt.tight_layout()

    if show:
        plt.show()

def plot_histograms(
    adata,
    obs_column: str,
    obs_values: list = None,
    var_names: list = None,
    bins: int = 30,
    alpha: float = 0.5,
    col_wrap: int = 3,
    figsize: tuple = (14, 6),
    show: bool = True
):
    """
    Plot histograms for selected genes stratified by a group in adata.obs.

    Parameters:
        adata (AnnData): AnnData object.
        obs_column (str): Column in adata.obs to group by.
        obs_values (list): Values in obs_column to include. If None, include all.
        var_names (list): Genes to plot. If None, use all.
        bins (int): Number of histogram bins.
        alpha (float): Opacity of histogram bars.
        col_wrap (int): Number of columns in FacetGrid layout.
        figsize (tuple): Figure size for each subplot row.
        show (bool): Whether to display the plot.

    Returns:
        None
    """
    if var_names is None:
        var_names = list(adata.var_names)
    
    if obs_values is not None:
        adata = adata[adata.obs[obs_column].isin(obs_values)].copy()

    obs_series = adata.obs[obs_column].astype(str)

    for var in var_names:
        if var not in adata.var_names:
            print(f"Warning: '{var}' not found in adata.var_names. Skipping.")
            continue

        df = pd.DataFrame({
            obs_column: obs_series.values,
            var: adata[:, var].X.toarray().flatten() if hasattr(adata[:, var].X, "toarray") else adata[:, var].X.flatten()
        })

        plt.figure(figsize=figsize)
        sns.histplot(
            data=df,
            x=var,
            hue=obs_column,
            bins=bins,
            kde=True,
            alpha=alpha,
            element="step",
            common_norm=False,
            palette="Set2"
        )
        plt.title(f"Histogram of {var} by {obs_column}")
        plt.tight_layout()
        if show:
            plt.show()
