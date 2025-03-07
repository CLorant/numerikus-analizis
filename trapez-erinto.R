LOWER_LIM <- 0
UPPER_LIM <- pi
N <- 10
X0 <- (LOWER_LIM + UPPER_LIM) / 2 

# Állítható függvény
f <- function(x) {
    x^3 * sin(x)
}

# Második derivált numerikus számítása
second_derivative <- function(f, x) {
    h <- 1e-5
    return((f(x + h) - 2 * f(x) + f(x - h)) / (h^2))
}

# Trapézszabály
trapezoidal_rule <- function(f, a, b, n) {
    x <- seq(a, b, length.out = n + 1)
    y <- f(x)
    h <- (b - a) / n
    integral <- h * (sum(y) - (y[1] + y[n + 1]) / 2)
    return(integral)
}

# Érintőformula (szakaszolt)
midpoint_rule <- function(f, a, b, n) {
    h <- (b - a) / n
    midpoints <- a - (h / 2) + (1:n) * h
    integral <- h * sum(f(midpoints))
    return(integral)
}

# Trapézszabály hibabecslés
trapezoidal_error_estimate <- function(f, a, b, n) {
    x_vals <- seq(a, b, length.out = 100)
    second_deriv_vals <- sapply(x_vals, function(x) second_derivative(f, x))
    max_second_deriv <- max(abs(second_deriv_vals))
    error_estimate <- ((b - a)^3 / (12 * n^2)) * max_second_deriv
    return(error_estimate)
}

# Trapézszabály utólagos hibabecslés
trapezoidal_a_posteriori_error <- function(f, a, b, n) {
    T_n <- trapezoidal_rule(f, a, b, n)
    T_2n <- trapezoidal_rule(f, a, b, 2 * n)
    error_estimate <- abs(T_n - T_2n)
    
    return(error_estimate)
}

# Érintőformula hibabecslés
midpoint_error_estimate <- function(f, a, b, n) {
    h <- (b - a) / n
    x_vals <- seq(a, b, length.out = 100)
    second_deriv_vals <- sapply(x_vals, function(x) second_derivative(f, x))
    max_second_deriv <- max(abs(second_deriv_vals))
    error_estimate <- ((b - a) * h^2 / 24) * max_second_deriv
    return(error_estimate)
}

get_ylim <- function(y_vals, pad_factor) {
    y_min <- min(y_vals)
    y_max <- max(y_vals)
    y_padding <- (y_max - y_min) * pad_factor
    ylim <- c(y_min, y_max + pad_factor)
}

plot_graph <- function(a, b, n, x0) {
    x_vals <- seq(a, b, length.out = 100)
    y_vals <- f(x_vals)
    ylim <- get_ylim(y_vals, 5)
    
    curve(f, a, b, col = "red", lwd = 2, main = "Trapézszabály és érintőformula", ylim=ylim)
    
    x_vals <- seq(a, b, length.out = n + 1)
    y_vals <- f(x_vals)
    for (i in 1:n) {
        polygon(c(x_vals[i], x_vals[i], x_vals[i+1], x_vals[i+1]),
                c(0, y_vals[i], y_vals[i+1], 0),
                col = rgb(0, 1, 0, alpha = 0.2), border = "green")
    }
    
    midpoint_x <- a + ((1:n) - 0.5) * (b - a) / n
    midpoint_y <- f(midpoint_x)
    points(midpoint_x, midpoint_y, col = "blue", pch = 16)
    lines(midpoint_x, midpoint_y, col = "blue", lwd = 2, lty = 2)
    
    trapezoidal_err <- trapezoidal_error_estimate(f, a, b, n)
    midpoint_err <- midpoint_error_estimate(f, a, b, n)
    a_posteriori_error <- trapezoidal_a_posteriori_error(f, a, b, n)
    
    legend_text <- c(
        paste("Trapézszabály hiba:", round(trapezoidal_err, 6)),
        paste("Érintőformula hiba:", round(midpoint_err, 6)),
        paste("Érintőformula utólagos hiba:", round(a_posteriori_error, 6))
    )
    
    y_left <- f(a)
    y_right <- f(b)
    
    if (y_left > y_right) {
        legend_pos <- "topright"
    }
    else {
        legend_pos <- "topleft"
    }
    
    legend(legend_pos, legend = c("f(x)", "Trapézszabály", "Érintőformula (Szakaszolt)", legend_text),
           col = c("red", "green", "blue", NA, NA, NA), 
           lwd = 2, lty = c(1, 1, 1, NA, NA, NA), 
           pch = c(NA, NA, 16, NA, NA, NA), 
           fill = c(NA, rgb(0, 0, 1, alpha = 0.2), NA, NA, NA, NA))
}

main <- function() {
    plot_graph(LOWER_LIM, UPPER_LIM, N, X0);
    
    numerical_integral <- trapezoidal_rule(f, LOWER_LIM, UPPER_LIM, N)
    midpoint_integral <- midpoint_rule(f, LOWER_LIM, UPPER_LIM, N)
    exact_integral <- integrate(f, LOWER_LIM, UPPER_LIM)$value
    
    trapezoidal_err <- trapezoidal_error_estimate(f, LOWER_LIM, UPPER_LIM, N)
    midpoint_err <- midpoint_error_estimate(f, LOWER_LIM, UPPER_LIM, N)
    a_posteriori_error <- trapezoidal_a_posteriori_error(f, LOWER_LIM, UPPER_LIM, N)
    percentage_trap_err = (trapezoidal_err / abs(exact_integral)) * 100
    percentage_midpoint_err = (midpoint_err / abs(exact_integral)) * 100
    percentage_a_posteriori_err = (a_posteriori_error / abs(exact_integral)) * 100
    
    cat("\n=== Paraméterek ===\n")
    cat("a =", LOWER_LIM, "\n")
    cat("b =", UPPER_LIM, "\n")
    cat("n =", N, "\n")
    cat("x0 =", X0, "\n")
    
    cat("\n=== Integrálási Eredmények ===\n")
    cat("Integrál értéke:\t", exact_integral, "\n")
    cat("Trapézszabály:\t\t", numerical_integral, "\n")
    cat("Érintőformula:\t\t", midpoint_integral, "\n")
    
    cat("\n=== Hibabecslések ===\n")
    cat("Trapézszabály:\t\t", trapezoidal_err, "\n")
    cat("Érintőformula:\t\t", midpoint_err, "\n")
    cat("Érintőformula utólagos:\t", a_posteriori_error, "\n")
    
    cat("\n=== Százalékos Eltérések ===\n")
    cat("Trapézszabály:\t\t", percentage_trap_err, "%\n")
    cat("Érintőformula:\t\t", percentage_midpoint_err, "%\n")
    cat("Érintőformula utólagos:\t", percentage_a_posteriori_err, "%\n")
}

main()
