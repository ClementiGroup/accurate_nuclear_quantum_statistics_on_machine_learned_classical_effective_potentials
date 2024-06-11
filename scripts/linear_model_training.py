import torch
import pickle
import numpy as np
from sklearn.model_selection import KFold
from sklearn.metrics import mean_squared_error
from utils import LinearRegressionForceMatching
from jsonargparse import CLI
from typing import Literal, Optional


def train_linear_fm(
    n_rbf: int,
    cutoff_lower: float,
    cutoff_upper: float,
    alpha: float,
    n_splits: int,
    model_type: str,  # Literal["linear_regression", "ridge"],
    splitting_random_state: int,
    position_file: str,
    force_file: str,
    records_to_skip: int = 0,
    device: str = "cpu",  # Literal['cuda', 'cpu'] = 'cpu',
    device_id: Optional[int] = None,
):
    """
    Train a linear model for force-matching

    Parameters:
    ----------
    n_rbf: int
        Number of radial basis functions
    cutoff_lower: float
        Lower cutoff
    cutoff_upper: float
        Upper cutoff
    alpha: float
        Regularization parameter of the ridge regression
    n_splits : int
        Number of splits in the dataset (number of folds, i.e, 5 for 5-fold cross-validation)
    model_type: str
        Type of the model passed to LinearRegressionForceMatching
    splitting_random_state: int
        Seed for splitting the dataset
    position_file : str
        Path to position file
    force_file : str
        Path to force file
    records_to_skip: int, default 0
        Number of records in the dataset to skip
    device: str, default 'cuda'
        device to use for model training
    device_id: int, optional
        device id
    """

    # Set the device
    device = torch.device(device, device_id)

    # Prepare training dataset
    input_data = np.load(position_file)[records_to_skip::]
    reference_data = np.load(force_file)[records_to_skip::]

    # Will do k-fold cross-validation.
    kf = KFold(n_splits=n_splits, random_state=splitting_random_state, shuffle=True)

    mse_training_list = []
    mse_test_list = []
    counter = 0  # Counter of folds
    for train_index, test_index in kf.split(input_data):
        print("Start new fold")
        X_train, X_test = input_data[train_index], input_data[test_index]
        y_train, y_test = reference_data[train_index], reference_data[test_index]
        training_input = torch.tensor(X_train, requires_grad=True)
        training_output = torch.tensor(y_train)

        test_input = torch.tensor(X_test, requires_grad=True)
        test_output = torch.tensor(y_test)

        model = LinearRegressionForceMatching(
            cutoff_lower=cutoff_lower,
            cutoff_upper=cutoff_upper,
            num_rbf=n_rbf,
            device=device,
            model_type=model_type,
            alpha=alpha,
        )
        model.to("cpu")
        model.fit(training_input, training_output)

        # Do prediction on two types of data:
        prediction_training = model.predict(training_input)

        mse_training = mean_squared_error(
            model.format_reference_forces(training_output), prediction_training
        )
        mse_training_list.append(mse_training)
        print("MSE training: ", mse_training)
        prediction = model.predict(test_input)
        mse = mean_squared_error(model.format_reference_forces(test_output), prediction)
        mse_test_list.append(mse)
        print("MSE test: ", mse)
        with open(f"optimized_model_fold_{counter}.pkl", "wb") as f:
            pickle.dump(model, f)
        counter += 1
    np.savetxt("optimal_params.txt", model.params)
    np.savetxt("training_mse.txt", np.array(mse_training_list))
    np.savetxt("test_mse.txt", np.array(mse_test_list))
    return


if __name__ == "__main__":
    CLI(train_linear_fm)
