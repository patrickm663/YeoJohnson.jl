# YeoJohnson

[![Build Status](https://github.com/patrickm663/YeoJohnson.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/patrickm663/YeoJohnson.jl/actions/workflows/CI.yml?query=branch%3Amain)

This package serves as a lightweight implementation of the [Yeo-Johnson](https://en.wikipedia.org/wiki/Power_transform#Yeo%E2%80%93Johnson_transformation) power transform.

## Algorithm
The Yeo-Johnson power transformation is an extension of the Box-Cox transformation that allows inputs that are not strictly positive. The intention behind applying a power transformation is to improve normality.

The formula for the transformation is as follows for each $y_{i}$ $\in$ $\mathbf{y}$:

$$\begin{equation}
\psi{(\lambda,y_{i})} =
    \begin{cases}
      ((y_{i}+1)^{\lambda}-1)/\lambda, & \text{if} \lambda \neq 0, y_i \ge 0 \\
      \log{(y_{i}+1)}, & \text{if } \lambda = 0, y_i \ge 0 \\
      -[(-y_{i}+1)^{2-\lambda}-1]/(2-\lambda), & \text{if}\ \lambda \neq 2, y_i < 0 \\
      -\log{(-y_{i}+1)}, & \text{if } \lambda = 2, y_i < 0 \\
    \end{cases}
    \end{equation}$$

The value $\lambda$ is estimated by minimising the negative likelihood function. The search interval is -2, 2 (by default).

## Installation
The project can be installed from the Julia REPL as follows once Package mode is active (pressing `]` in the REPL):
```julia
add "https://github.com/patrickm663/YeoJohnson.jl"
```
The project's chief dependency is the [Optim.jl](https://github.com/JuliaNLSolvers/Optim.jl) package which estimates for $\lambda$.

## Function
The only exported function currently is `yeojohnson()`.

It accepts the following parameters:
- `x`: a vector to be transformed
- `min`: the minimum value in the search range for lambda (DEFAULT = -2)
- `max`: the maximum value in the search range for lambda (DEFAULT = 2)
- `\lambda`: a user-defined lambda value if the optimiser is not being used
- `opt`: a boolean value for whether to oprtimise for lambda or use the user-defined value

## Example

```julia
using YeoJohnson
using Distribution, Random # To create dummy data

d = Gamma(3, 5)
x = rand(d, 10_000)

y = yeojohnson(x)
```
![Untransformed Gamma](/images/untransformed_gamma.png)

![Transformed Gamma](/images/transformed_gamma.png)

## Source
The original source of the algorithm is Yeo, In-Kwon and Johnson, Richard A. (2000). "A New Family of Power Transformations to Improve Normality or Symmetry". Biometrika. 87 (4): 954–959. doi:10.1093/biomet/87.4.954

## License
This project is licensed under the MIT License.
