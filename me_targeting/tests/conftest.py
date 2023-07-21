# Copyright 2017 BioInformatics Research Center, KAIST
import pytest
from os.path import join, abspath, dirname

data_model_dir = join(dirname(abspath(__file__)), 'data')

@pytest.fixture(scope="function")
def model_file():
    model_file = join(data_model_dir, 'iAF1260_COBRA.xml')
    return model_file