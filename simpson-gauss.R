library(pracma)

f <- function(x) {
    x^3 * sin(x)
}

simpson_simple <- function(f, a, b) {
    h <- (b - a) / 2
    integral <- (h / 3) * (f(a) + 4 * f((a + b) / 2) + f(b))
    return(integral)
}

simpson_composite <- function(f, a, b, n) {
    if (n %% 2 != 0) {
        stop("Az n értékének párosnak kell lennie.")
    }
    
    h <- (b - a) / n
    x <- seq(a, b, length.out = n + 1)
    y <- f(x)
    
    integral <- (h / 3) * (y[1] + 4 * sum(y[seq(2, n, by = 2)]) + 2 * sum(y[seq(3, n-1, by = 2)]) + y[n + 1])
    return(integral)
}

fourth_derivative <- function(f, x) {
    h <- 1e-5
    deriv <- (f(x + 2*h) - 4*f(x + h) + 6*f(x) - 4*f(x - h) + f(x - 2*h)) / (h^4)
    if (!is.finite(deriv)) {
        deriv <- 0
    }
    return(deriv)
}

simpson_simple_error_estimate <- function(f, a, b) {
    h <- (b - a) / 2
    x_vals <- seq(a, b, length.out = 100)
    fourth_deriv_vals <- sapply(x_vals, function(x) fourth_derivative(f, x))
    max_fourth_deriv <- max(abs(fourth_deriv_vals), na.rm = TRUE)
    error_estimate <- (max_fourth_deriv * (b - a)^5) / 2880
    return(error_estimate)
}

simpson_error_estimate <- function(f, a, b, n) {
    h <- (b - a) / n
    x_vals <- seq(a, b, length.out = 100)
    fourth_deriv_vals <- sapply(x_vals, function(x) fourth_derivative(f, x))
    max_fourth_deriv <- max(abs(fourth_deriv_vals), na.rm = TRUE)
    error_estimate <- (max_fourth_deriv * (b - a)^5) / (2880 * n^4)
    return(error_estimate)
}

simpson_a_posteriori_error <- function(f, a, b, n) {
    S_n <- simpson_composite(f, a, b, n)
    S_2n <- simpson_composite(f, a, b, 2 * n)
    error_estimate <- abs(S_n - S_2n)
    return(error_estimate)
}

gauss_quadrature <- function(f, a, b, n) {
    nodes_weights <- gaussLegendre(n, a, b)
    nodes <- nodes_weights$x
    weights <- nodes_weights$w
    integral <- sum(weights * f(nodes))
    return(integral)
}

gauss_error_estimate <- function(f, a, b, n) {
    x_vals <- seq(a, b, length.out = 100)
    deriv_vals <- sapply(x_vals, function(x) {
        h <- 1e-5
        (f(x + 2*h) - 4*f(x + h) + 6*f(x) - 4*f(x - h) + f(x - 2*h)) / (h^4)
    })
    max_deriv <- max(abs(deriv_vals), na.rm = TRUE)
    error_estimate <- (max_deriv * (b - a)^(2*n + 3)) / (factorial(2*n + 2))
    return(error_estimate)
}

plot_graph <- function(a, b, n, f, simpson_composite_result, a_posteriori_error, gauss_result, gauss_err) {
    curve(f, from = a, to = b, col = "black", lwd = 2, main = "Simpson-formula és Gauss-kvadratúra")
    grid(nx = NULL, ny = NULL, lty = 2, col = "lightgray", lwd = 2)
    
    simpson_x <- seq(a, b, length.out = n + 1)
    simpson_y <- f(simpson_x)
    lines(simpson_x, simpson_y, col = "blue", lwd = 2, lty = 2)
    points(seq(a, b, length.out = n + 1), f(seq(a, b, length.out = n + 1)), col = "blue", pch = 16)
    
    nodes_weights <- gaussLegendre(n, a, b)
    gauss_x <- nodes_weights$x
    gauss_y <- f(gauss_x)
    lines(gauss_x, gauss_y, col = "red", lwd = 2, lty = 3)
    points(gauss_x, gauss_y, col = "red", pch = 17)
    
    mid_x <- (a + b) / 2
    y_left <- f(a)
    y_right <- f(b)
    
    if (y_left > y_right) {
        legend_pos <- "topright"
    }
    else {
        legend_pos <- "topleft"
    }
    
    legend(legend_pos, legend = c(
        "f(x)",
        "Simpson-formula",
        "Gauss-kvadratúra",
        paste("Simpson utólagos hiba:", round(a_posteriori_error, 6)),
        paste("Gauss-kvadratúra hiba:", round(gauss_err, 6))
        
    ), col = c("black", "blue", "red", NA, NA), lwd = 2, lty = c(1, 2, 3, NA, NA), pch = c(NA, 16, 17, NA, NA))
}

main <- function() {
    a <- 0
    b <- pi
    n <- 10
    
    simpson_simple_result <- simpson_simple(f, a, b)
    simpson_composite_result <- simpson_composite(f, a, b, n)
    exact_integral <- integrate(f, a, b)$value
    
    a_posteriori_error <- simpson_a_posteriori_error(f, a, b, n)
    percentage_a_posteriori_err <- (a_posteriori_error / abs(exact_integral)) * 100
    
    gauss_result <- gauss_quadrature(f, a, b, n)
    gauss_err <- gauss_error_estimate(f, a, b, n)
    percentage_gauss_err <- (gauss_err / abs(exact_integral)) * 100
    
    cat("\n=== Simpson-formula ===\n")
    cat("Integrál értéke:\t\t", exact_integral, "\n")
    cat("Simpson-formula (egyszerű):\t", simpson_simple_result, "\n")
    cat("Simpson-formula (összetett):\t", simpson_composite_result, "\n")
    cat("Utólagos hiba:\t\t\t", a_posteriori_error, "\n")
    cat("Százalékos hiba (utólagos):\t", percentage_a_posteriori_err, "%\n")
    
    cat("\n=== Gauss-kvadratúra ===\n")
    cat("Gauss-kvadratúra:\t\t", gauss_result, "\n")
    cat("Hiba (analitikus):\t\t", gauss_err, "\n")
    cat("Százalékos hiba (analitikus):\t", percentage_gauss_err, "%\n")
    
    plot_graph(a, b, n, f, simpson_composite_result, a_posteriori_error, gauss_result, gauss_err)
}

main()