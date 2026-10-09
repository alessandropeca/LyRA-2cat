
def LYRA_sigma_part1(args):
	globals().update(args)
	global save_output
	global save_images
	global path_output
	global filename_outinfo
	global r_in_tm
	global r_min_nm
	global r_max_nm
	global r_lr
	global ra_input
	global dec_input
	global ra_output
	global dec_output
	global sigma_input
	global sigma_output
	global input_cat_type
	global output_cat_type
	global ID_input
	global ID_output
	global mags_input
	global input_cat_type
	global noxy
	global nomag
	global sigma_input
	global path_images
	global save_images
	global needed_only_once
	global area_input
	global area_output

	if needed_only_once == False:
		import numpy as np
		import LYRA_functions as lyra

		print('\nLoading errors...')
		#-------------------------------------------------------------------------------------- input:
		if input_cat_type == 0:
			sigma_input = sigma_input
			print('Sigma input:', round(sigma_input,2))
		elif input_cat_type == 1:
			sigma_input = lyra.quadratic_sum(ra_input_err[ra_input_err!=noxy],dec_input_err[dec_input_err!=noxy])
			print('Sigma input: using single positional errors')
		elif input_cat_type == 2:
			sigma_input = lyra.quadratic_sum(ra_input_err, dec_input_err)
			print('Sigma mean input:', round(sigma_input,2))

		#-------------------------------------------------------------------------------------- output:
		if output_cat_type == 0:
			sigma_output = sigma_output
			print('Sigma output:', round(sigma_output,2))
		elif output_cat_type == 1:
			sigma_output = lyra.quadratic_sum(ra_output_err[ra_output_err!=noxy],dec_output_err[dec_output_err!=noxy])
			print('Sigma output: using single positional errors.')
		elif output_cat_type == 2:
			sigma_output = lyra.quadratic_sum(ra_output_err, dec_output_err)
			print('Sigma mean output:', round(sigma_output,2))


		print('... done.\n')
		return globals()
