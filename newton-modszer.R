# Deriv csomag betöltése
library(Deriv)

# Függvény a y értékek tartományának meghatározásához
get_range <- function(y_vals) {
    list(y_min = min(y_vals), y_max = max(y_vals))
}

# Függvény a grafikon előkészítéséhez
plot_prepare <- function(name, f, f_prime, x0, epsilon = 1e-6, max_iter = 100) {
    root <- x0

    # Intervallum pontosítás a newton-módszer prototípus futtatásával
    for (i in 1:max_iter) {
        root <- root - f(root) / f_prime(root)
        if (abs(f(root)) < epsilon) break
    }
    
    zoom_factor <- abs(f(root))
    if (zoom_factor < epsilon) {
        zoom_factor <- 1
    }
    
    a <- root - zoom_factor
    b <- root + zoom_factor
    
    x_vals <- seq(a, b, length.out = 100)
    y_vals <- sapply(x_vals, f)
    range <- get_range(y_vals)
    
    # Plot the function
    plot(x_vals, y_vals, type = "l", col = "blue", lwd = 2, xlab = "x", ylab = "y", 
         ylim = c(range$y_min, range$y_max), main = paste("Newton-módszer |", name))
    abline(h = 0, col = "black")  # Add a horizontal line at y = 0
    
    return(range$y_min)
}

# Függvény az iterációs lépések ábrázolásához
plot_iter <- function(iter, f, x, y, y_min) {
    lines(c(x, y), c(f(x), 0), col = "red", lwd = 2)

    if (iter == 1) {
        lines(c(x, x), c(f(x), y_min), col = "black", lwd = 1, lty = 2)
    }

    # TODO: Change from fixed number iter check to showing the dotted line and the x label if they don't overlap
    if (iter < 3) {
        lines(c(y, y), c(f(y), y_min), col = "black", lwd = 1, lty = 2)
    }
    if (iter < 4) {
        text(x, y_min, labels = paste0("x", iter - 1), pos = 1, col = "red", cex = 0.8)
    }
}

# Függvény a végső ábrázoláshoz
plot_finalize <- function(y, y_min) {
    points(y, 0, col = "darkblue", pch = 18, lwd = 3, cex = 2.5)

    legend("topleft", legend = c("y = 0", "y = f(x)", "Iterációk", "Gyök"), 
            col = c("black", "blue", "red", "darkblue"), 
            lty = c(1, 1, 1, NA), lwd = c(1, 2, 1, NA), pch = c(NA, NA, NA, 18), 
            pt.cex = c(NA, NA, NA, 2.5))
}

# Függvény az iterációs lépések kiírásához
print_step <- function(iter, x, y) {
    cat("\nIteráció", iter, "\n")
    cat("x:", x, "\n")
    cat("y:", y, "\n")
}

# Newton-módszer algoritmus
newton_method <- function(name, f, x0, epsilon = 1e-6, max_iter = 100) {
    f_prime <- Deriv(f)

    y_min <- plot_prepare(name, f, f_prime, x0, epsilon, max_iter)
    
    x <- x0
    x_prev <- x0
  
    for (iter in 1:max_iter) {
        y <- x - f(x) / f_prime(x)
        
        plot_iter(iter, f, x, y, y_min)
        print_step(iter, x, y)

        if (abs(y - x) <= epsilon && abs(f(y)) <= epsilon) {
            plot_finalize(y, y_min)
            
            return(list(root = y, iterations = iter))
        }
        
        x_prev <- x
        x <- y
    }
  
    cat("\nMaximum iteráció elérve. A gyök nem talált.\n")
    return(list(root = x, iterations = max_iter))
}

read_choice <- function() {
    message("1: Lineáris konvergencia");
    message("2: Négyzetes konvergencia");
    message("3: Oszcillációs konvergencia");
    message("4: Lassú konvergencia");
    
    choice <- 1
    
    repeat {
        choice <- as.numeric(readline())
        
        if (!is.na(choice) && choice %in% 1:4) {
            return(choice)
        }
        
        cat("Érvénytelen érték, próbáld újra!\n")
    }
}

get_demo <- function(choice) {
    demos <- list(
        "1" = list(name = "Lineáris konvergencia",
                    f = function(x){ (x + 2) / 3 }, x0 = 0),
        "2" = list(name = "Négyzetes konvergencia",
                    f = function(x){ x^2 - 2 }, x0 = 1),
        "3" = list(name = "Oszcillációs konvergencia",
                    f = function(x){ cos(x) }, x0 = 1),
        "4" = list(name = "Lassú konvergencia",
                    f = function(x){ x^10 - 1 }, x0 = 0.5)
    )
	
    demo <- demos[[as.character(choice)]]
}

get_function_body <- function(f) {
    gsub("\\s|[{}]", "", deparse(body(f)))
}

main <- function() {
    repeat {
        choice <- read_choice()
        demo = get_demo(choice)

		cat("\n=== ",  demo$name ," ===\n")
        cat("f(x) =", get_function_body(demo$f), "\n")
		cat("x0 = ", demo$x0, "\n")

        result <- newton_method(demo$name, demo$f, demo$x0)
        cat("\nGyök: ", result$root, "\n")
        cat("Iterációk száma: ", result$iterations, "\n")

        message("\nÚjra futtatod? i / n")
        if (tolower(readline()) != "i") break
    }
}

main()
