## 1.0.1.2

Added an optional ridge term `lambda`. If `lambda > 0`, an identity matrix scaled by `lambda` is added to the covariance matrices of all components both during initialization and estimation to ensure positive definiteness.
