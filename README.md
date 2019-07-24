Rama GMM
--------

Code to fit a basic Gaussian Mixture Model to Ramachandran plots.

![Trained GMM](images/fullclassifications-1.png)

## Code organization

1. `ramachandran_training.Rmd`: Train the models (slow). Outputs tsv files.
2. `ramachandran_model.Rmd`: Apply the models. Main report
3. `ramachandran_helpers.R`: helper functions

## Output

1. `model_full.tsv` Parameters for the trained model
2. `ramachandran_model.pdf` Report with full details

## Training Data

Phi-psi angles from good quality protein structures were used to train the model. These were downloaded from the [Protein Geometry Database](http://pgd.cgrb.oregonstate.edu/). For details see `ramachandran_model.pdf`

## License

Copyright Spencer Bliven

Licensed under the Unlicense (public domain).

