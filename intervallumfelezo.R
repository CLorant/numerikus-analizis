f <- function(x) {
    x^2 - 50
}

get_ymin <- function(a, b) {
    x_values <- seq(a, b, length.out = 1)
    y_values <- f(x_values)
    y_min <- min(y_values)
    
    return(y_min)
}

get_ymax <- function(a, b) {
    x_values <- seq(a, b, length.out = 1)
    y_values <- f(x_values)
    y_max <- max(y_values)
    
    return(y_max)
}

plot_prepare <- function(a, b, y_min, y_max) {
    plot(f, from=a-1, to=b+1, col="blue", lwd=2, main="Intervallumfelező módszer",
         xlab="x-tengely", ylab="y-tengely", ylim=c(y_min, y_max))
    abline(h=0, col="black", lwd=2)
}

plot_step <- function(a, b, y_min, segment_height, max_iter) {
    # Függőleges vonalak
    y_vertical_a <- seq(y_min, segment_height, length.out = 2)
    x_vertical_a <- rep(a, length(y_vertical_a))
    lines(x_vertical_a, y_vertical_a, lwd=2, col="gray")
    
    y_vertical_b <- seq(y_min, segment_height, length.out = 2)
    x_vertical_b <- rep(b, length(y_vertical_b))
    lines(x_vertical_b, y_vertical_b, lwd=2, col="gray")
    
    # Vízszintes vonal
    x_horizontal <- seq(a, b, length.out = 2)
    y_horizontal <- rep(segment_height - 2, length(x_horizontal))
    lines(x_horizontal, y_horizontal, col="red", lwd=2)
    
    if (max_iter < 8) {
        text((a + b) / 2, y_horizontal + 40 / max_iter, labels = paste0("[", a, ", ", b, "]"), col="black", cex=0.8)
    }
}

print_step <- function(iter, a, b, c, fa, fc, max_iter) {
    cat("\nIteráció", iter, "\n")
    cat("a:\t", a, "\n")
    cat("b:\t", b, "\n")
    cat("c:\t", c, "\n")
    cat("f(a):\t", fa, "\n")
    cat("f(c):\t", fc, "\n")
    
    if (iter == max_iter) {
        cat("\nEredmények:");
    }
    else if (fa * fc < 0) {
        cat("Folytatás a(z) [", a, ", ", c, "] intervallumon ([a, c])\n")
    }
    else {
        cat("Folytatás a(z) [", c, ", ", b, "] intervallumon ([c, b])\n")
    }
}

# Intervallumfelező algoritmus
bisection_method <- function(f, a, b, tol = 1e-10, max_iter = 5) {
    segment_height = (max_iter + 1) * 50
    y_min <- get_ymin(a, b)
    y_max <- get_ymax(a, b) + segment_height + 100
    plot_prepare(a, b, y_min, y_max)
    
    iter <- 1
    fa <- f(a)
    
    while ((b - a) > 2 * tol && iter <= max_iter) {
        plot_step(a, b, y_min, segment_height, max_iter)
        segment_height <- segment_height - 50
        
        c <- (a + b) / 2
        fc <- f(c)
        
        print_step(iter, a, b, c, fa, fc, max_iter)
        
        if (fa * fc < 0) {
            b <- c
        } else {
            a <- c
            fa <- fc
        }
        
        iter <- iter + 1
    }

    # Utolsó lépés kiírása
    plot_step(a, b, y_min, segment_height, max_iter)

    legend("topleft", legend = c("f(x)", "Vizsgált intervallum", "Gyök"),
           col = c("blue", "red", "black"), 
           lwd = 2, lty = c(1, 1, NA), 
           pch = c(NA, NA, 18))

    root = (a + b) / 2
    points(root, 0, pch=18, lwd=3, cex=2.5)

    return(root)
}

main <- function() {
    a <- 5
    b <- 8.5
    
    if (f(a) * f(b) >= 0) {
        cat("Hiba: Az intervallum nem tartalmaz gyököt, vagy több gyök is lehet benne.\n")
        a <- 0
        b <- 10
        cat("Az intervallum kibővítve erre: [", a, ", ", b, "]\n")
        
        if (f(a) * f(b) >= 0) {
            stop("Hiba: Az új intervallum sem tartalmaz gyököt.")
        }
    }
    
    if (a > b) {
        cat("Figyelem: a > b, ezért felcseréljük a két végpontot.\n")
        temp <- a
        a <- b
        b <- temp
    }
    
    approx_root <- bisection_method(f, a, b)
    exact_root <- uniroot(f, interval = c(a, b), tol = 1e-10)$root
    diff <- abs(exact_root - approx_root)
    diff_percent <- (diff / abs(exact_root)) * 100
    
    cat("\nKözelítő gyök:\t\t", approx_root, "\n")
    cat("Pontos gyök:\t\t", exact_root, "\n")
    cat("Abszolút eltérés:\t", diff, "\n")
    cat("Százalékos eltérés:\t", diff_percent, "%\n")
}

main()
