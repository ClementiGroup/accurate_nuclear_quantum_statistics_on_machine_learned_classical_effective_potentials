root=/group/ag_clementi_cmb/users/iryna/PI_CG/bulk_h2o/mace_model2/trial_1
EPOCH=15


for folder in $root/temp_300*; do
	    # Get the last 4 symbols of the folder name
		seed=${folder: -4}
		# get the model 
		for model in $folder/ckpt/epoch\=${EPOCH}-validation_loss=*; do
			# get the basename of the model
			model_name=$(basename $model)
			# append seed to the model name
			model_name=seed_${seed}_${model_name}
			# copy the model to the models folder
			cp $model ./300/$model_name
		done
	done
