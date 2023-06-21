import pandas as pd
import subprocess
import pytest

def tblite_xtbml_subprocess(dir,bin,args):
    stat = subprocess.call(
    [bin] +["run"]+["coord.xyz"]+ ["--xtbml_xyz"]+args,
    shell=False,
    stdin=None,
    stderr=subprocess.DEVNULL,cwd=dir
)

import os
def test_dir(dir,ref_file,bin,args):
    
    df_ref = pd.read_csv(ref_file)
    dir = os.path.dirname(ref_file)
    if "xyz" in dir:
        to_be_dropped = ["w_dens","w_alpha"]
    else:
        #in the old xyz print out there was error for extended Z only and e only print out
        to_be_dropped = ["w_dens","w_alpha","delta_qm_Z_xx","delta_qm_Z_yy","delta_qm_Z_xy","delta_qm_e_xx","delta_qm_e_yy","delta_qm_e_xy"]
    df_ref = df_ref.drop(to_be_dropped,axis=1)
    
    cwd = os.getcwd()
    bin = cwd+"/"+bin
    tblite_xtbml_subprocess(dir,bin,args)
    df_test = pd.read_csv(dir+"/ml_feature_tblite.csv")

    to_be_dropped_test = ['delta_qm_e_xx', 'q_d', 'qm_s_zz', 'chem.pot', 'qm_s_yz', 
    'delta_qm_Z_yz', 'dipm_p', 'qm_s_xz', 'qm_s_yy', 'dipm_s_z', 'delta_dipm_Z', 'dipm_d', 
    'qm_s_xy', 'qm_d', 'delta_qm_e', 'delta_qm_e_xz', 'qm_d_zz', 'qm_d_xx',
    'dipm_d_z', 'dipm_d_y', 'qm_p_yz', 'qm_p',"delta_qm_e_yz",
    'qm_p_xz', 'dipm_A', 'dipm_s_x', 'dipm_p_x', 'qm_p_zz', 'response', 'qm_p_xy',
    'q_p', 'qm_s', 'dipm_s_y', 'gap', 'qm_d_xy', 'dipm_d_x', 'qm_s_xx', 'delta_qm_e_zz',
    'delta_dipm_A', 'qm_d_yz','qm_A', 'qm_d_yy', 'dipm_p_y', 'q_s', 'dipm_p_z', 'delta_qm_Z', 'qm_p_xx', 
    'delta_qm_Z_zz', 'dipm_s', 'delta_qm', 'qm_d_xz', 'delta_qm_Z_xz', 'qm_p_yy', 'delta_dipm_e']
    df_test = df_test.drop(to_be_dropped_test,axis=1)
    df = abs(df_ref-df_test )
    thr = 0.0002
    filter_ = (df > thr).any(axis=0)
    print(df_test["HOAO_a"],df_ref["HOAO_a"])
    print(df.loc[:,filter_])
    os.remove(dir+"/ml_feature_tblite.csv")
    
    assert df.loc[:,filter_].empty
        

