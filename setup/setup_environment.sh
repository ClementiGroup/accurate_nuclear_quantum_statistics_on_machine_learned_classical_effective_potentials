conda env create --file environment.yml 

source ~/.bashrc
conda activate picg
pip install   torch==1.12.0+cu113 -f https://download.pytorch.org/whl/torch_stable.html
