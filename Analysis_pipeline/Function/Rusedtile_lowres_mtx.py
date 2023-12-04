
'''
This script is for ST data SME modified
'''
#!/usr/bin/env python
# conding:utf-8 -*-
from pathlib import Path, PurePath
from typing import Union, Dict, Optional, Tuple, BinaryIO

import argparse
import stlearn as st
import scanpy as sc
import numpy as np
from numpy import random,mat
from pathlib import Path
import pandas as pd
from scipy import io,sparse
import os

import h5py
import json
from matplotlib.image import imread
import anndata
from anndata import (
    AnnData,
    read_csv,
    read_text,
    read_excel,
    read_mtx,
    read_loom,
    read_hdf,
)
from anndata import read as read_h5ad
import logging as logg


def read_visium_lowres_mtx(
    path: Union[str, Path],
    genome: Optional[str] = None,
    *,
    count_file: str = "filtered_feature_bc_matrix/",
    library_id: str = None,
    load_images: Optional[bool] = True,
    source_image_path: Optional[Union[str, Path]] = None,
) -> AnnData:
    path = Path(path)
    adata = sc.read_10x_mtx(path / count_file)

    adata.uns["spatial"] = dict()

    # from h5py import File
    # 
    # with File(path / count_file, mode="r") as f:
    #     attrs = dict(f.attrs)
    # if library_id is None:
    #     library_id = str(attrs.pop("library_ids")[0], "utf-8")

    adata.uns["spatial"][library_id] = dict()

    if load_images:
        files = dict(
            tissue_positions_file=path / 'spatial/tissue_positions_list.csv',
            scalefactors_json_file=path / 'spatial/scalefactors_json.json',
            hires_image=path / 'spatial/tissue_hires_image.png',
            lowres_image=path / 'spatial/tissue_lowres_image.png',
        )

        # check if files exists, continue if images are missing
        for f in files.values():
            if not f.exists():
                if any(x in str(f) for x in ["hires_image", "lowres_image"]):
                    logg.warning(
                        f"You seem to be missing an image file.\n"
                        f"Could not find '{f}'."
                    )
                else:
                    raise OSError(f"Could not find '{f}'")

        adata.uns["spatial"][library_id]['images'] = dict()
        for res in ['hires', 'lowres']:
            try:
                adata.uns["spatial"][library_id]['images'][res] = imread(
                    str(files[f'{res}_image'])
                )
            except Exception:
                raise OSError(f"Could not find '{res}_image'")

        # read json scalefactors
        adata.uns["spatial"][library_id]['scalefactors'] = json.loads(
            files['scalefactors_json_file'].read_bytes()
        )

        # adata.uns["spatial"][library_id]["metadata"] = {
        #     k: (str(attrs[k], "utf-8") if isinstance(attrs[k], bytes) else attrs[k])
        #     for k in ("chemistry_description", "software_version")
        #     if k in attrs
        # }

        # read coordinates
        positions = pd.read_csv(files['tissue_positions_file'], header=None)
        positions.columns = [
            'barcode',
            'in_tissue',
            'array_row',
            'array_col',
            'pxl_col_in_fullres',
            'pxl_row_in_fullres',
        ]
        positions.index = positions['barcode']

        adata.obs = adata.obs.join(positions, how="left")

        adata.obsm['spatial'] = adata.obs[
            ['pxl_row_in_fullres', 'pxl_col_in_fullres']
        ].to_numpy()
        adata.obs.drop(
            columns=['barcode', 'pxl_row_in_fullres', 'pxl_col_in_fullres'],
            inplace=True,
        )

        # put image path in uns
        if source_image_path is not None:
            # get an absolute path
            source_image_path = str(Path(source_image_path).resolve())
            adata.uns["spatial"][library_id]["metadata"]["source_image_path"] = str(
                source_image_path
            )

    return adata
  
def Read10X_lowres_mtx(
    path: Union[str, Path],
    genome: Optional[str] = None,
    count_file: str = "filtered_feature_bc_matrix",
    library_id: str = None,
    load_images: Optional[bool] = True,
    quality= "lowres",
    image_path: Union[str, Path] = None,
) -> AnnData:

    adata = read_visium_lowres_mtx(
        path,
        genome=genome,
        count_file=count_file,
        library_id=library_id,
        load_images=load_images,
    )
    adata.var_names_make_unique()

    if library_id is None:
        library_id = list(adata.uns["spatial"].keys())[0]

    if quality == "fulres":
        image_coor = adata.obsm["spatial"]
        img = plt.imread(image_path, 0)
        adata.uns["spatial"][library_id]["images"]["fulres"] = img
    else:
        scale = adata.uns["spatial"][library_id]["scalefactors"][
            "tissue_" + quality + "_scalef"
        ]
        image_coor = adata.obsm["spatial"] * scale

    adata.obs["imagecol"] = image_coor[:, 0]
    adata.obs["imagerow"] = image_coor[:, 1]
    adata.uns["spatial"][library_id]["use_quality"] = quality

    return adata
  
def ME_normalize(inDir,outDir,sample):
    print (sample, "start SME normalize")

    #read data
    data=Read10X_lowres_mtx(path = inDir)
    data.var_names_make_unique()
    data.layers['raw_count']=data.X
    #tile data
    TILE_PATH=Path(os.path.join(outDir,'{0}_tile'.format(sample)))
    TILE_PATH.mkdir(parents=True,exist_ok=True)
    
    #tile morphology
    st.pp.tiling(data,TILE_PATH,crop_size=40)
    st.pp.extract_feature(data)

    ###process data
    st.pp.normalize_total(data)
    st.pp.log1p(data)

    #gene pca dimention reduction
    st.em.run_pca(data,n_comps=50,random_state=0)
    
    #stSME to normalise log transformed data
    st.spatial.SME.SME_normalize(data, use_data="raw",weights =  "weights_matrix_gd_md")

    #convert SME_norm data to sparesmatrix
    raw_SME_normalized = mat(data.obsm['raw_SME_normalized'])
    raw_SME_normalizedA = sparse.csr_matrix(raw_SME_normalized)
    print ("matrix convert ok!")
    
    io.mmwrite(os.path.join(outDir,'{0}_raw_SME_normalizeA.mtx'.format(sample)),raw_SME_normalizedA)
    print("Morphology adjusted ok!")

    return raw_SME_normalizedA
   
