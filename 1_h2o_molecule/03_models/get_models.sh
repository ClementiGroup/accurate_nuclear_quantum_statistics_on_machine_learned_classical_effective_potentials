root=/import/a12/users/iryna/PI_CG/h2o_molecule/mace_model
EPOCH=60

for trial in {1..6}; do
	echo $trial
	for folder in $root/trial_$trial/temp_${trial}00*; do
	    # Get the last 4 symbols of the folder name
		seed=${folder: -4}
		# get the model 
		for model in $folder/ckpt/epoch\=${EPOCH}-validation_loss=*; do
			# get the basename of the model
			model_name=$(basename $model)
			# append seed to the model name
			model_name=seed_${seed}_${model_name}
			# copy the model to the models folder
			cp $model ./${trial}00/$model_name
		done
	done
done
