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
    s_n <- simpson_composite(f, a, b, n)
    s_2n <- simpson_composite(f, a, b, 2 * n)
    error_estimate <- abs(s_n - s_2n)
    return(error_estimate)
}

legendre_poly <- function(n, x) {
    if (n == 0) {
        return(rep(1, length(x)))
    }
    else if (n == 1) {
        return(x)
    }
    else {
        p0 <- rep(1, length(x))
        p1 <- x
        for (k in 2:n) {
            pk <- ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
            p0 <- p1
            p1 <- pk
        }
        return(p1)
    }
}

legendre_poly_deriv <- function(n, x) {
    if (n == 0) {
        return(rep(0, length(x)))
    }
    else {
        return(n * (x * legendre_poly(n, x) - legendre_poly(n - 1, x)) / (x^2 - 1))
    }
}

find_legendre_roots <- function(n, a, b, tol = 1e-10, max_iter = 100) {
    roots <- numeric(n)
    for (i in 1:n) {
        x0 <- cos(pi * (i - 0.25) / (n + 0.5))
        for (iter in 1:max_iter) {
            Pn <- legendre_poly(n, x0)
            Pn_deriv <- legendre_poly_deriv(n, x0)
            x1 <- x0 - Pn / Pn_deriv
            if (abs(x1 - x0) < tol) {
                break
            }
            x0 <- x1
        }
        roots[i] <- x1
    }
    return(roots)
}

gauss_weights <- function(n, a, b, roots) {
    weights <- numeric(n)
    for (i in 1:n) {
        xi <- roots[i]
        Pn_deriv <- legendre_poly_deriv(n, xi)
        weights[i] <- 2 / ((1 - xi^2) * Pn_deriv^2)
    }
    return(weights)
}

gauss_quadrature <- function(f, a, b, n) {
    roots <- find_legendre_roots(n, a, b)
    weights <- gauss_weights(n, a, b, roots)
    
    nodes <- ((b - a) * roots + (a + b)) / 2
    integral <- sum(weights * f(nodes)) * (b - a) / 2
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

get_ylim <- function(y_vals, pad_factor) {
    y_min <- min(y_vals)
    y_max <- max(y_vals)
    y_padding <- (y_max - y_min) * pad_factor
    ylim <- c(y_min, y_max + pad_factor)
}

plot_graph <- function(a, b, n, f, simpson_composite_result, a_posteriori_error, gauss_result, gauss_err) {
    x_vals <- seq(a, b, length.out = 100)
    y_vals <- f(x_vals)
    ylim <- get_ylim(y_vals, 5)
    
    curve(f, from = a, to = b, col = "black", lwd = 2, main = "Simpson-formula és Gauss-kvadratúra (Legendre polinommal)", ylim=ylim)
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
        
    ),
    col = c("black", "blue", "red", NA, NA),
    lwd = 2, lty = c(1, 2, 3, NA, NA),
    pch = c(NA, 16, 17, NA, NA))
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