# MatrixNormalShrinkFlop
Fast shrinkage based estimator for matrix normal covariance

## Usage 

FlipFlop.R contains the algorithm for the estimator, which is called regsinkhorn. Based on the forthcoming joint work https://arxiv.org/abs/2110.07583 with Akshay Ramachandran, Rafael Oliveira, and Michael Walter and http://proceedings.mlr.press/v125/franks20a/franks20a.pdf with Ankur Moitra.

The file ShrinkFlop.rmd compares this algorithm with other estimators such as KGlasso and Gemini in various generative models.
