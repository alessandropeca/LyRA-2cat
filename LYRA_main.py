import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import LYRA_functions as lyra
import time
import argparse
import math
import ast
import sys
import warnings
import IPython
warnings.filterwarnings("ignore")

parser = argparse.ArgumentParser("thx xPecax")
# NOTE this is when we send in the catalogs
# NOTE angular resolution only comes from the secondary (worse angluar resolution, i.e., highest)
parser.add_argument("ang_res", type=float)
# positional errors comes from both primary and secondary
parser.add_argument("pos_err_primary")
parser.add_argument("pos_err_secondary")
parser.add_argument("--mean_err", default=None)
parser.add_argument("file_primary", help='primary catalog')
parser.add_argument("file_secondary", help='secondary catalog')
# they all get changed to the same value
# parser.add_argument("--r_in_tm", type=float, default=2)
# parser.add_argument("--r_lr", type=float, default=5)
# parser.add_argument("--r_min_nm", type=int, default=5)
# # NOTE experimental: for now we'll do 5 * minimum ang resolution
# parser.add_argument("--r_max_nm", type=int, default=30)
# leave it
parser.add_argument("--distrib_bins", type=int, default=21)
parser.add_argument("--path_LR", default="./") # main table/result
parser.add_argument("--path_output", default = None) # txt file with info on what done by LYRA
parser.add_argument("--path_images", default = None) 
# NOTE we'll convert our tables to -99 to be compatible with this
parser.add_argument("--noxy", type=int, default=-99)
parser.add_argument("--nomag", type=int, default=-99)
parser.add_argument("--save_images", default=False)
parser.add_argument("--save_output", default=True) # txt file with info on what done by LYRA
# this is the most important output info
parser.add_argument("--filename_LR", default="tmp_table") # main table/result
# gives some metadata as well
parser.add_argument("--filename_outinfo", default="tmp_info") # txt file with info on what done by LYRA
# magnitude column comes from the input
parser.add_argument("--mag_col", default=None)
# NEW parameters
parser.add_argument("--min_sample_size", default=4) # minimum input_cat sample size to start LYRA (4: silverman et al 1986 'required sample size for given accuracy')
parser.add_argument("--thresh_perc_mag", default=0.1) # minimum % of good mag values in a given column to build the mag prior.
parser.add_argument("--bayes", default=True) # Use bayesian statistics
parser.add_argument("--add_prior", default=1) # if you want to add another baysian factor (--bayes True is required here)
parser.add_argument("--compl_prior", default=True) # Q_int value
parser.add_argument("--ID_prim", default="ID") # take also different colum names
parser.add_argument("--RA_prim", default="RA")
parser.add_argument("--DEC_prim", default="DEC")
parser.add_argument("--ID_sec", default="ID")
parser.add_argument("--RA_sec", default="RA")
parser.add_argument("--DEC_sec", default="DEC")
parser.add_argument("--rel_threshold", default=0.2) # reliability threshold for R_flag
parser.add_argument("--convhull_area", default=True) # less realibile but speed up the process
# parser.add_argument("--source_names", default="auto") # "auto" tries to read src names from catalogs, otherwise it uses indeces; "method" uses an external file to assign names
parser.add_argument("--custom_RC_threshold", default=-99) # if you don't want to find the maximum in RC but check for possible plateaux
parser.add_argument("--add_name_output", default="")
toc = time.time()
print()
goto_pos_match = False
exit_table_only = False
mag_failure = False

args = parser.parse_args()
if args.path_output is None:
    args.path_output = args.path_LR + "output/"
if args.path_images is None:
    args.path_images = args.path_LR + "output/images/"

globals().update(args.__dict__)

# for some reason we need this, it seems that parsing boolean values does not work. So:
if (bayes == "False") or (bayes == False):
    bayes = False
else:
    bayes = True
if (save_images == "False") or (save_images == False):
    save_images = False
else:
    save_images=True
if (save_output == "False") or (save_output == False):
    save_output = False
else:
    save_output = True
if (compl_prior == "False") or (compl_prior == False):
    compl_prior = False
else:
    compl_prior = True
if (convhull_area == "False") or (convhull_area == False):
    convhull_area = False
else:
    convhull_area = True

if (custom_RC_threshold != -99) or (custom_RC_threshold != '-99'):
    custom_RC_threshold = float(custom_RC_threshold)

if (nomag != -99) or (nomag != '-99'):
    nomag = float(nomag)

if save_output == True:
    outinfo = open(path_output+filename_outinfo,"w")
    outinfo.write('###################\n'+'### OUTPUT INFO ###\n' + '###################\n\n')
    outinfo.close()
    outinfo.close()


if save_output == True:
    parser.add_argument("ang_res", type=float)
    outinfo = open(path_output+filename_outinfo,"a")
    outinfo.write("\nRequired parameters: ")
    outinfo.write("\nfile_primary: " + str(file_primary))
    outinfo.write("\nfile_secondary: " + str(file_secondary))
    outinfo.write("\nang_res: " + str(ang_res))
    outinfo.write("\npos_err_primary: " + str(pos_err_primary))
    outinfo.write("\npos_err_secondary: " + str(pos_err_secondary))
    outinfo.write("\nOptional parameters: ")
    outinfo.write("\nsave_images: " + str(save_images))
    outinfo.write("\nsave_output: " + str(save_output))
    outinfo.write("\nbayes: " + str(bayes))
    outinfo.write("\nmean_err: " + str(mean_err))
    outinfo.write("\ndistrib_bins: " + str(distrib_bins))
    outinfo.write("\npath_LR: " + str(path_LR))
    outinfo.write("\npath_output: " + str(path_output))
    outinfo.write("\npath_images: " + str(path_images))
    outinfo.write("\nnoxy: " + str(noxy))
    outinfo.write("\nnomag: " + str(nomag))
    outinfo.write("\nfilename_LR: " + str(filename_LR))
    outinfo.write("\nfilename_outinfo: " + str(filename_outinfo))
    outinfo.write("\nmag_col: " + str(mag_col))
    outinfo.write("\nmin_sample_size: " + str(min_sample_size))
    outinfo.write("\nthresh_perc_mag: " + str(thresh_perc_mag))
    outinfo.write("\nadd_prior: " + str(add_prior))
    outinfo.write("\ncompl_prior: " + str(compl_prior))
    outinfo.write("\nID_prim: " + str(ID_prim))
    outinfo.write("\nRA_prim: " + str(RA_prim))
    outinfo.write("\nDEC_prim: " + str(DEC_prim))
    outinfo.write("\nID_sec: " + str(ID_sec))
    outinfo.write("\nRA_sec: " + str(RA_sec))
    outinfo.write("\nDEC_sec: " + str(DEC_sec))
    outinfo.write("\nrel_threshold: " + str(rel_threshold))
    outinfo.write("\nconvhull_area: " + str(convhull_area))
    outinfo.write("\ncustom_RC_threshold: " + str(custom_RC_threshold))
    # outinfo.write("\nsource_names: " + str(source_names))
    outinfo.close()

if lyra.is_number(pos_err_primary) == True:
    input_cat_type = 0
    pos_err_primary = float(pos_err_primary)
else:
    if mean_err != None:
        input_cat_type = 2
    elif lyra.is_number(pos_err_primary) == False:
        input_cat_type = 1
    else:
        print("\n###### Error: Cannot read the input_cat_type/output_cat_type. ######")
        sys.exit()

if lyra.is_number(pos_err_secondary) == True:
    output_cat_type = 0
    pos_err_secondary = float(pos_err_secondary)
else:
    if mean_err != None:
        output_cat_type = 2
    elif lyra.is_number(pos_err_secondary) == False:
        output_cat_type = 1
    else:
        print("\n###### Error: Cannot read the input_cat_type/output_cat_type. ######")
        sys.exit()

r_in_tm = ang_res
r_lr = ang_res
r_min_nm = ang_res
# NOTE experimental: for now we'll do 5 * minimum ang resolution
r_max_nm = 10 * r_min_nm

# to read CSV but also FITS/FIT tables
if (file_primary.lower().endswith(('.fits')) == True) or (file_primary.lower().endswith(('.fit')) == True):
    from astropy.table import Table
    dat_in = Table.read(file_primary, format='fits')
    primary_df = dat_in.to_pandas()
else:
    primary_df = pd.read_csv(file_primary)



##### cut for mag columns if all of them are -99 -> for S82-XL only #####
if mag_col is not None:
    mag_cols_ = mag_col.split(",")
    mag_cols_[-1] = mag_cols_[-1].strip() # needed in some linux systems
    primary_df = primary_df[~(primary_df[mag_cols_] == nomag).all(axis=1)]
    primary_df.reset_index(drop=True, inplace=True)
#####



if (file_secondary.lower().endswith(('.fits')) == True) or (file_secondary.lower().endswith(('.fit')) == True):
    from astropy.table import Table
    dat_out = Table.read(file_secondary, format='fits')
    secondary_df = dat_out.to_pandas()
else:
    secondary_df = pd.read_csv(file_secondary)


# Check for possible duplicates (RA, DEC):
# print("\nprimary cols:")
plist = primary_df.columns.to_list()
# print(plist)
# IPython.embed()
if primary_df.duplicated(subset=(RA_prim,DEC_prim)).any() == True:
    # primary_df.drop_duplicates(subset=(RA_prim,DEC_prim), inplace=True, ignore_index=True)
    # only for AGN-DB!!! otherwise use the line above
    # agg_dict = {col: 'first' for col in primary_df.columns}
    
    # substring = "mag"
    # prim_col_list = []
    # for ls in plist:
    #     if substring in ls:
    #         prim_col_list.append(ls)
    #     if ls == 'pos_error_Temporary':
    #         prim_col_list.append(ls)
    # agg_dict = {col: 'mean' for col in prim_col_list}
    # agg_dict['ID_new'] = 'first'
    # agg_dict['ID'] = 'first'
    # agg_dict['matched_attr_tuples'] = 'first'
    # print(agg_dict)
    # primary_df = primary_df.groupby([RA_prim,DEC_prim]).agg(agg_dict).reset_index()
    primary_df.drop_duplicates(subset=(RA_prim, DEC_prim), inplace=True, ignore_index=True)
    ####################################################

    print("\n### Warning: Found duplicates in primary catalog (based on RA and DEC) "+str(file_primary)+". Now merged. ###")
    if save_output == True:
        lyra.print_output(path_output+filename_outinfo, "### Warning: Found duplicates in primary catalog (based on RA and DEC) "+str(file_primary)+". Now merged. ###")

ra_input = np.array(primary_df[RA_prim])
dec_input = np.array(primary_df[DEC_prim])

# if source_names == "auto":
# 	try:
# 		ID_input = np.array(primary_df[ID_prim])
# 	except:
# 		ID_input = np.array(primary_df.index)
# 		print("Primary catalog: indeces used as src names.")
# 		if save_output == True:
# 			lyra.print_output(path_output+filename_outinfo, "Primary catalog: indeces used as src names.")
# elif source_names == "method":
# 	try:
# 		ID_input = lyra.LYRA_source_names(ra_input, dec_input)
# 	except:
# 		print("\n### Error: Cannot read source names. ###")
# 		if save_output == True:
# 			lyra.print_output(path_output+filename_outinfo, "### Error: Cannot read source names. ###")
# else:
# 	print("\n### Error: Cannot read source names. ###")
# 	if save_output == True:
# 		lyra.print_output(path_output+filename_outinfo, "### Error: Cannot read source names. ###")

# ID_input contains indeces used by LYRA for the match, while ID_input_names are the source names if present
try:
    # ID_input_names = np.array(primary_df[ID_prim]) # commented for AB-DB only!!!!!!!
    ID_input = np.array(primary_df.index)
    ID_input_names = ID_input # this is for AGN-DB only!!!!!!!
except:
    ID_input = np.array(primary_df.index)
    ID_input_names = ID_input
    ID_prim = "ID_prim"
    primary_df[ID_prim] = ID_input
    print("Primary catalog: indeces used as src names.")
    if save_output == True:
        lyra.print_output(path_output+filename_outinfo, "Primary catalog: indeces used as src names.")


# Check for possible duplicates (RA, DEC):
# print("secondary cols:")
slist = secondary_df.columns.to_list()
# print(slist)
if secondary_df.duplicated(subset=(RA_sec,DEC_sec)).any() == True:
    # secondary_df.drop_duplicates(subset=(RA_sec,DEC_sec), inplace=True, ignore_index=True)
    # only for AGN-DB!!! otherwise use the line above
    # substring = "mag"
    # sec_col_list = []
    # for ls in slist:
    #     if substring in ls:
    #         sec_col_list.append(ls)
    #     if ls == 'pos_error_Temporary':
    #         sec_col_list.append(ls)
    agg_dict = {col: 'first' for col in secondary_df.columns 
               if col not in [RA_sec, DEC_sec]}
    # agg_dict = {col: 'first' for col in sec_col_list}
    # agg_dict['ID_new'] = 'first'
    # agg_dict['ID'] = 'first'
    # agg_dict['matched_attr_tuples'] = 'first'
    # print(agg_dict)
    secondary_df = secondary_df.groupby([RA_sec,DEC_sec]).agg(agg_dict).reset_index()
    ####################################################
    print("\n### Warning: Found duplicates in secondary catalog (based on RA and DEC) "+str(file_primary)+". Now dropped. ###")
    if save_output == True:
        lyra.print_output(path_output+filename_outinfo, "### Warning: Found duplicates in secondary catalog (based on RA and DEC) "+str(file_primary)+". Now dropped. ###")

ra_output = np.array(secondary_df[RA_sec])
dec_output = np.array(secondary_df[DEC_sec])

try:
    # ID_output_names = np.array(secondary_df[ID_sec]) # commented for AB-DB only!!!!!!!
    ID_output = np.array(secondary_df.index)
    ID_output_names = ID_output # this is for AB-DB only!!!!!!!
except:
    ID_output = np.array(secondary_df.index)
    ID_output_names = ID_output
    ID_sec = "ID_sec"
    secondary_df[ID_sec] = ID_output
    print("Secondary catalog: indeces used as src names.")
    if save_output == True:
        lyra.print_output(path_output+filename_outinfo, "Secondary catalog: indeces used as src names.")


if "origin" in primary_df.columns:
    pri_origin = primary_df["origin"].tolist()
    print("\n\n HERE\n pri_origin", pri_origin)
    print(pri_origin)
    print(primary_df["origin"])
    if not isinstance(primary_df["origin"][0], tuple):
        primary_df['origin'] = primary_df['origin'].apply(ast.literal_eval)
        if not isinstance(primary_df["origin"][0], tuple):
            print("The value is not a tuple.")
    pri_origin = primary_df["origin"].tolist()
else:
    primary_df["origin"] = [(1,) for _ in range(len(primary_df))]
    pri_origin = primary_df["origin"].tolist()

if "origin" in secondary_df.columns:
    sec_origin = secondary_df["origin"].tolist()
    print("\n\n HERE\n sec_origin", sec_origin)
    if not isinstance(secondary_df["origin"][0], tuple):
        secondary_df['origin'] = secondary_df['origin'].apply(ast.literal_eval)
        if not isinstance(secondary_df["origin"][0], tuple):
            print("The value is not a tuple.")
    sec_origin = secondary_df["origin"].tolist()
else:
    secondary_df["origin"] = [(2,) for _ in range(len(secondary_df))]
    sec_origin = secondary_df["origin"].tolist()

# # # check if it is not the first time matching  
# if "origin" in primary_df.columns:
#     pri_origin = primary_df["origin"].tolist()
# else:
#     file_pri_tw_a = file_primary.split("/")[-1]
#     if file_pri_tw_a[0] == 'd':
#         file_pri_tw_ = file_pri_tw_a.split(".csv")[0]
#         file_pri_tw = int(file_pri_tw_.split("d")[-1])
#         pri_origin = [(file_pri_tw,) for _ in range(len(primary_df))]
#     else:
#         file_pri_tw_b = file_pri_tw_a.split("lyra_")[-1]
#         file_pri_tw = int(file_pri_tw_b.split(".csv")[0])
#         pri_origin = [(file_pri_tw,) for _ in range(len(primary_df))]
# print("pri_origin", pri_origin)

# for i, element in enumerate(pri_origin):
#     if not isinstance(element, tuple):
#         print(f"The element at index {i} is not a tuple.")
#         print(f"Type: {type(element)}, Value: {element}")
# assert all(isinstance(element, tuple) for element in pri_origin), "Not all elements are tuples."

# if "origin" in secondary_df.columns:
#     sec_origin = secondary_df["origin"].tolist()
# else:
#     file_sec_tw_a = file_secondary.split("/")[-1]
#     if file_sec_tw_a[0] == 'd':
#         file_sec_tw_ = file_sec_tw_a.split(".csv")[0]
#         file_sec_tw = int(file_sec_tw_.split("d")[-1])
#         sec_origin = [(file_sec_tw,) for _ in range(len(secondary_df))]
#     else:
#         file_sec_tw_b = file_sec_tw_a.split("lyra_")[-1]
#         file_sec_tw = int(file_sec_tw_b.split(".csv")[0])
#         sec_origin = [(file_sec_tw,) for _ in range(len(secondary_df))]
# print("sec_origin", sec_origin)

# for i, element in enumerate(sec_origin):
#     if not isinstance(element, tuple):
#         print(f"The element at index {i} is not a tuple.")
#         print(f"Type: {type(element)}, Value: {element}")
# assert all(isinstance(element, tuple) for element in sec_origin), "Not all elements are tuples."

# input/output catalogs positional errors: 0 = static number for all sources, 1 = single error for each source, 2 = mean error for all sources
if input_cat_type == 0:
    sigma_input = pos_err_primary
elif input_cat_type == 1:
    ra_input_err = np.array(primary_df[pos_err_primary]) / np.sqrt(2)
    dec_input_err = np.array(primary_df[pos_err_primary]) / np.sqrt(2)
elif input_cat_type == 2:
    err_in_tmp = np.array(primary_df[pos_err_primary])
    ra_input_err = np.mean(err_in_tmp[err_in_tmp>0] / np.sqrt(2))
    dec_input_err = np.mean(err_in_tmp[err_in_tmp>0] / np.sqrt(2))

if output_cat_type == 0:
    sigma_output = pos_err_secondary
    ra_output_err = pos_err_secondary / np.sqrt(2)
    dec_output_err = pos_err_secondary / np.sqrt(2)
elif output_cat_type == 1:
    ra_output_err = np.array(secondary_df[pos_err_secondary]) / np.sqrt(2)
    dec_output_err = np.array(secondary_df[pos_err_secondary]) / np.sqrt(2)
elif output_cat_type == 2:
    err_out_tmp = np.array(secondary_df[pos_err_secondary])
    ra_output_err = np.mean(err_out_tmp[err_out_tmp>0] / np.sqrt(2))
    dec_output_err = np.mean(err_out_tmp[err_out_tmp>0] / np.sqrt(2))

# create folders
if not os.path.isdir(path_output):
    os.mkdir(path_output)
if save_images == True:
    if not os.path.isdir(path_images):
        os.mkdir(path_images)

if (len(np.unique(ra_input[ra_input != noxy])) > min_sample_size) and (len(np.unique(dec_input[dec_input != noxy])) > min_sample_size):
    if (len(np.unique(ra_output[ra_output != noxy])) > min_sample_size) and (len(np.unique(dec_output[dec_output != noxy])) > min_sample_size):

        # if more than one mag_in columuns, then iterate the code on all the mag_in columns
        all_lr_mag_in = []
        all_fr_r = []
        all_bfactor = []
        all_prior = []
        all_lr_LR = []
        all_LR_threshold = []
        all_rel = []
        all_F_flag = []

        needed_only_once = False
        if mag_col is not None:
            eff_mag_cols = []
            mag_cols = mag_col.split(",")
            mag_cols[-1] = mag_cols[-1].strip() # needed in some linux systems

            for smag in mag_cols:
                mag_input = np.array(primary_df[smag])
                if np.median(mag_input[(mag_input!= nomag)]) < 1: # we have fluxes! try with -2.5log10(flux) ~ magnitudes
                    mag_input = -2.5*np.log10(mag_input)

                #--- check magnitudes and coords ---
                mag_input = np.nan_to_num(mag_input) # check for nans values
                mag_input[(mag_input<5) | (mag_input>40)] = nomag # check for weird values
                # mag_input[(mag_input>20)] = nomag

                if (len(mag_input[(mag_input!= nomag)]) >= min_sample_size) and ((len(mag_input[mag_input>0])/len(mag_input) >= thresh_perc_mag) or (len(mag_input[mag_input>0]) >= 1000)): # we need enough good entries
                    from LYRA_LR_compl_part4 import LYRA_LR_compl_part4
                    globals().update(LYRA_LR_compl_part4(globals()))

                    if goto_pos_match == False:
                        all_lr_mag_in = np.concatenate((all_lr_mag_in, lr_mag_in))
                        all_bfactor = np.concatenate((all_bfactor, bfactor))
                        all_prior = np.concatenate((all_prior, prior))
                        all_lr_LR = np.concatenate((all_lr_LR, lr_LR))
                        all_rel = np.concatenate((all_rel, rel))
                        eff_mag_cols.append(smag)
                        all_LR_threshold = np.append(all_LR_threshold, LR_threshold)
                
                else: # check for enough good values
                    print("\n### Warning: Not enough good values for mag_input: " + smag + ", so I skip it. ###")
                    if save_output == True:
                        lyra.print_output(path_output+filename_outinfo, "### Warning: Not enough good values for mag_input: " + smag + ", so I skip it. ###")
            
            if len(all_lr_mag_in) > 0:
                all_lr_mag_in = all_lr_mag_in.reshape(len(eff_mag_cols), len(lr_mag_in))
                all_bfactor = all_bfactor.reshape(len(eff_mag_cols), len(lr_mag_in))
                all_prior = all_prior.reshape(len(eff_mag_cols), len(lr_mag_in))
                all_lr_LR = all_lr_LR.reshape(len(eff_mag_cols), len(lr_mag_in))
                all_rel = all_rel.reshape(len(eff_mag_cols), len(lr_mag_in))

            else:
                goto_pos_match = True
                print("\n### Warning: No mag_input columns with enough good values. ###\n")
                if save_output == True:
                    lyra.print_output(path_output+filename_outinfo, "### Warning: No mag_input columns with enough good values. ###")
        else:
            mag_cols = []
            goto_pos_match = True
            print("\n### Warning: No mag_input given, go to simple positional match. ###\n")

        if goto_pos_match == False: # we have magnitude columns at this point
            
            data_to_write = {"OutID": outID_names,
                            # "OutID": lr_outID,
                            "RAout": lr_outRA,
                            "DECout": lr_outDEC,
                            "InID": inID_names,
                            # "InID": lr_inID,
                            "RAin": lr_inRA,
                            "DECin": lr_inDEC,
                            "r": lr_r,
                            "match_err": sigma_match,
                            "origin": origin,
                            "S_flag": 0,
                            "C_flag": 0
                            }
            
            for m in range(len(eff_mag_cols)):
                data_to_write[eff_mag_cols[m]] = all_lr_mag_in[m]
                data_to_write['Bfact_'+eff_mag_cols[m]] = all_bfactor[m]
                data_to_write['Prior_'+eff_mag_cols[m]] = all_prior[m]
                data_to_write['Post_'+eff_mag_cols[m]] = all_lr_LR[m]
                data_to_write['Prel_'+eff_mag_cols[m]] = all_rel[m]
                data_to_write['F_flag_'+eff_mag_cols[m]] = -99

            # F_flag (Flavour_flag): 1 = best match with LR>=LRth, 2 = secondary match with LR>=LRth and rel>=rel_threshold, 0 = other cases LR<LRth, -99 = no match.
            # C_flag: 1 = best match, 0 = low prob. or no match.
            print('\nComputing flags...')
            df = pd.DataFrame(data_to_write)
            df = df.sort_values(by=["OutID"])
            df_ok = df[df["InID"]!= -99]
            s_outID = np.unique(df_ok['OutID'])

            if len(eff_mag_cols) > 1:
                index=0
                for subid in s_outID:
                    index+=1
                    cc_flag = False
                    df_sub = df[df['OutID']==subid]
                    df_sub_indeces = df_sub.index.to_list()
                    
                    for m in range(len(eff_mag_cols)):
                        lr_tmp = df_sub['Post_'+eff_mag_cols[m]]
                        prel_tmp = df_sub['Prel_'+eff_mag_cols[m]]
                        local_LRth = all_LR_threshold[m]
                        
                        # first, compute F_flag for each eff_mag_cols[m]
                        if len(df_sub[df_sub['Post_'+eff_mag_cols[m]] >= local_LRth]) > 0: # at least one Post value above the threshold
                            cc_flag = True
                            idx_list = df_sub.index.to_list()
                            for idx in idx_list:
                                max_lr = max(lr_tmp)
                                if max_lr > 0:
                                    if (lr_tmp[idx] == max_lr) and (lr_tmp[idx] >= local_LRth):
                                        df['F_flag_'+eff_mag_cols[m]][idx] = 1
                                        df_sub['F_flag_'+eff_mag_cols[m]][idx] = 1
                                    else:
                                        if (lr_tmp[idx] >= local_LRth) and (prel_tmp[idx] >= rel_threshold):
                                            df['F_flag_'+eff_mag_cols[m]][idx] = 2
                                            df_sub['F_flag_'+eff_mag_cols[m]][idx] = 2
                                        else: # pos match only
                                            df['F_flag_'+eff_mag_cols[m]][idx] = 0 
                                            df_sub['F_flag_'+eff_mag_cols[m]][idx] = 0
                        
                        else: # for df_subs with only pos match with Post<LR
                            if len(df_sub[df_sub['Post_'+eff_mag_cols[m]] > 0]) > 0: # at least one != -99
                                idx_list = df_sub.index.to_list()
                                for idx in idx_list:
                                    if df_sub['Post_'+eff_mag_cols[m]][idx] > 0:
                                        df['F_flag_'+eff_mag_cols[m]][idx] = 0 
                                        df_sub['F_flag_'+eff_mag_cols[m]][idx] = 0

                    # second, compute C_flag: from here, only one inID source in df['OutID']==subid has to exit with C_flag = 1
                    if cc_flag == True:
                        save_idx = [] # indeces with at least a possible final C_flag == 1
                        save_count = []
                        
                        for idx_m in df_sub_indeces:
                            idx_line_f = []
                            for m in range(len(eff_mag_cols)):
                                idx_line_f.append(df_sub['F_flag_'+eff_mag_cols[m]][idx_m]) # same line, m F_flags
                                # S_flag: flag for possible secondary matches. ==1 if at least one time there is F_flag = 1 (that can happen multiple times)
                                if df_sub['F_flag_'+eff_mag_cols[m]][idx_m] == 1:
                                    df["S_flag"][idx_m] = 1

                            fcount = idx_line_f.count(1) # how many F_flag=1
                            save_count.append(fcount)
                            save_idx.append(idx_m)
                        
                        mm=0
                        for m in range(len(save_count)):
                            if save_count[m] == max(save_count):
                                mm+=1
                        
                        if mm == 1: # only one best line with more 1s, when comparing to the other lines
                            for m in range(len(save_count)):
                                if save_count[m] == max(save_count):
                                    if df["S_flag"][save_idx[m]] != 0:
                                        df['C_flag'][save_idx[m]] = 1
                        elif mm > 1: # more than one line with same max number of 1s -> check on reliabilities
                            possible_idx = []
                            for m in range(len(save_count)):
                                if save_count[m] == max(save_count):
                                    possible_idx.append(save_idx[m]) # check reliabilities on these indeces
                        
                            lines_max_rel = []
                            for p_idx in possible_idx:
                                line_rel = []
                                for m in range(len(eff_mag_cols)):
                                    if df_sub['F_flag_'+eff_mag_cols[m]][p_idx] == 1:
                                        line_rel.append(df_sub['Prel_'+eff_mag_cols[m]][p_idx])
                                    else:
                                        line_rel.append(0)
                                
                                local_max_rel = max(line_rel) # local (same idx) maxim rel
                                try:
                                    if len(local_max_rel) > 0: # more than one Rel value with the same max value -> pick one
                                        local_max_rel = local_max_rel[0]
                                except:
                                    local_max_rel = local_max_rel
                                lines_max_rel.append(local_max_rel) # same line max rel
                            
                            max_rel = max(lines_max_rel) # multiple lines (same outID) max rel -> C_flag = 1
                            for mm in range(len(lines_max_rel)): # line_max_rel and possible_idx have the same lenght
                                if lines_max_rel[mm] == max_rel: # if not unique, same value -> pick one
                                    best_idx = possible_idx[mm]
                            
                            if df["S_flag"][best_idx] != 0:
                                df['C_flag'][best_idx] = 1 # otherwise it is already = 0
                        
                    prog = 100*index/len(s_outID)
                    print("[%-20s] %d%%" % ('='*(int(prog/5)), (int(prog))), end="\r")

            elif len(eff_mag_cols) == 1:
                index=0
                for subid in s_outID:
                    index+=1
                    df_sub = df[df['OutID']==subid] # df_sub has same indeces of df
                    lr_tmp = df_sub['Post_'+eff_mag_cols[0]]
                    prel_tmp = df_sub['Prel_'+eff_mag_cols[0]]
                    local_LRth = all_LR_threshold[0]

                    if len(df_sub[df_sub['Post_'+eff_mag_cols[0]] >= local_LRth]) > 0:
                        idx_list = df_sub.index.to_list()
                        for idx in idx_list:
                            max_lr = max(lr_tmp)
                            if (lr_tmp[idx] == max_lr) and (lr_tmp[idx] >= local_LRth):
                                df['F_flag_'+eff_mag_cols[0]][idx] = 1
                                df['S_flag'][idx] = 1
                                df['C_flag'][idx] = 1
                            else:
                                if (lr_tmp[idx] >= local_LRth) and (prel_tmp[idx] >= local_LRth):
                                    df['F_flag_'+eff_mag_cols[0]][idx] = 2
                                else:
                                    df['F_flag_'+eff_mag_cols[0]][idx] = 0
                                df['C_flag'][idx] = 0
                
                    prog = 100*index/len(s_outID)
                    print("[%-20s] %d%%" % ('='*(int(prog/5)), (int(prog))), end="\r")
            
            else:
                print("\n Something went wrong: len(eff_mag_cols) <= 0. Go to positional match.")
                if save_output == True:
                    lyra.print_output(path_output+filename_outinfo, "Something went wrong: len(eff_mag_cols) <= 0. Go to positional match.")
                goto_pos_match = True
            
            
            df = df.sort_values(by=["OutID"])
            print('\n... done.\n')

            print('Cleaning multiple in ID counterpars...')
            df_sc = df.sort_values(by=["OutID"])
            df_sc_ok = df_sc[df_sc["InID"]!= -99]
            s_inID = np.unique(df_sc_ok['InID'])

            check_multiple_inID = df_sc_ok[df_sc_ok['C_flag'] == 1]['InID'].duplicated().any() # returns True if more than one inID with C_flag == 1
            if check_multiple_inID == True:
                df_tmp = df_sc[df_sc['C_flag']==1]
                df_sc_ok = df_sc[df_sc["InID"]!= -99]
                ss_inID = np.unique(df_tmp['InID'])
                index=0
                for subid in s_inID:
                    index+=1
                    df_sub = df_sc[df_sc['InID']==subid]
                    df_sub = df_sub[df_sub['C_flag']==1]
                    cf_tmp = list(df_sub['C_flag'])
                    b_index = df_sub['C_flag'].index.tolist() # indeces of same inID with C_flag=1
                    
                    if cf_tmp.count(1) > 1: # if same inID is associated with more than one ouID and C_flag=1
                        prel_max = []
                        b_index2 = []
                        for mags_name in eff_mag_cols:
                            b_index2.append(df_sub[df_sub['Prel_'+mags_name] == max(df_sub['Prel_'+mags_name])].index)
                            prel_max.append(max(df_sub['Prel_'+mags_name])) # max(Prel) of each mag col
                            
                        best_max = max(prel_max) # global rel max
                        best_index = b_index2[prel_max.index(max(prel_max))] # global rel max index
                        try:
                            if len(best_index) > 1: # it means we have two entries with the same max value
                                best_index = best_index[0] # so both the sources are ok, take one
                        except:
                            best_index = best_index
                        
                        for idx in range(len(b_index)):
                            if (b_index[idx] == best_index):
                                outid_tmp = df_sc['OutID'][b_index[idx]]
                                df_sub2 = df_sc[(df_sc['OutID'] == outid_tmp) & (df_sc['C_flag'] == 1)]
                                if len(df_sub2) <= 1:
                                    df_sc['C_flag'][b_index[idx]] = 1
                                else:
                                    df_sc['C_flag'][b_index[idx]] = 0
                            else: # these are the indeces that have C_flag=1 and need to be C_flag=0
                                df_sc['C_flag'][b_index[idx]] = 0
                                outid_tmp = df_sc['OutID'][b_index[idx]]
                                change_idx = False
                                # check if there is some source that can be turned into C_flag=0
                                df_sub2 = df_sc[(df_sc['OutID'] == outid_tmp)]
                                df_sub2_idx = df_sc[(df_sc['OutID'] == outid_tmp)].index
                                tmp_f_idx = []
                                tmp_p = []
                                for ind in df_sub2_idx:
                                    if ind != b_index[idx]:
                                        if (df_sub2['InID'][ind] in ss_inID) == False: # the inID must not already be taken
                                            # check for F_flag=1 or F_flag=2
                                            for mags_name in eff_mag_cols:
                                                if (df_sub2['F_flag_'+mags_name][ind] == 1) or (df_sub2['F_flag_'+mags_name][ind] == 2):
                                                    tmp_f_idx.append(ind)
                                                    tmp_p.append(df_sub2['Prel_'+mags_name][ind])
                                                    change_idx = True
                                                    
                                for mr in range(len(tmp_f_idx)):
                                    if tmp_p[mr] == max(tmp_p): # if more than one with the same max, pick one
                                        idx_to_reinclude = tmp_f_idx[mr]
                                if change_idx == True:
                                    df_sc['C_flag'][idx_to_reinclude] = 1
                                    ss_inID = np.append(ss_inID, df_sc['InID'][idx_to_reinclude])

                    prog = 100*index/len(s_inID)
                    print("[%-20s] %d%%" % ('='*(int(prog/5)), (int(prog))), end="\r")
            
            
            print('\n... done.')
            
            # order columns
            column_to_move = df_sc.pop('S_flag')
            df_sc["S_flag"] = column_to_move
            column_to_move = df_sc.pop('C_flag')
            df_sc["C_flag"] = column_to_move

            # set origin of the matched catalog
            assert all(isinstance(element, tuple) for element in df_sc["origin"].to_list()), "Not all elements are tuples."
            # df_sc["origin"] = [file_pri_tw + file_sec_tw if inid != -99 else file_sec_tw for inid in df_sc["InID"]]
            # assert all(isinstance(element, tuple) for element in df_sc["origin"].to_list()), "Not all elements are tuples."

            df = df_sc.sort_values(by=["OutID"])
            best_matches = df['C_flag'].tolist()
            print("\n\nPositional matches:", d)
            print("Best matches:", best_matches.count(1))
            if save_output == True:
                lyra.print_output(path_output+filename_outinfo, "-----------------------------------\nPositional matches: " + str(d))
                lyra.print_output(path_output+filename_outinfo, "Best matches: " + str(best_matches.count(1)))

            
            # write LYRA outputs
            print('\nWriting outputs...')
            out_file_name = path_output+filename_LR+".csv"
            df.to_csv(out_file_name, sep=',', index=False, header=True)
            print('... done.\n')

    else:
        goto_pos_match = True
        bayes = False
        needed_only_once = False
        print("\n### Warning: Output entries are not enough for building likelihood probabilities. ###")
        print("bayes is now:", bayes)
        if save_output == True:
            lyra.print_output(path_output+filename_outinfo, "### Warning: Output entries are not enough for building likelihood probabilities. ###")
            lyra.print_output(path_output+filename_outinfo, "bayes is now: " + str(bayes) + "\n")

else:
    goto_pos_match = True
    bayes = False
    needed_only_once = False
    print("\n### Warning: Input entries are not enough for building likelihood probabilities. ###")
    print("bayes is now:", bayes)
    if save_output == True:
        lyra.print_output(path_output+filename_outinfo, "### Warning: Output entries are not enough for building likelihood probabilities. ###")
        lyra.print_output(path_output+filename_outinfo,"bayes is now: " + str(bayes) + "\n")

if goto_pos_match == True:
    # call a single positional matching
    if bayes == True:
        from LYRA_pos_match_only_bayes import LYRA_pos_match_only_bayes
        LYRA_pos_match_only_bayes(globals())
    else:
        from LYRA_pos_match_only import LYRA_pos_match_only
        LYRA_pos_match_only(globals())

# time calculation
tic = time.time()
tt = tic - toc
tt_str = 'seconds.'
if tt > 60:
    tt = tt/60
    tt_str = 'minutes.'
    if tt > 3600:
        tt = tt/3600
        tt_str = 'hours.'
print('Computational time:', round(tt,1), tt_str, '\n')
if save_output == True:
    lyra.print_output(path_output+filename_outinfo, "Computational time: " + str(round(tt,1)) + " " + tt_str)

print('\nCheck duplicates [ID] in '+path_output+filename_LR+'.csv :')
out_dupl, in_dupl = lyra.check_duplicates_ids(path_output+filename_LR+".csv")
print('Out duplicates:', out_dupl, "from", file_secondary)
print('In duplicates:', in_dupl, "from", file_primary)
print('\nCheck duplicates [RA,DEC] in '+path_output+filename_LR+'.csv :')
out_dupl, in_dupl = lyra.check_duplicates_radec(path_output+filename_LR+".csv")
print('Out duplicates:', out_dupl, "from", file_secondary)
print('In duplicates:', in_dupl, "from", file_primary)
print()