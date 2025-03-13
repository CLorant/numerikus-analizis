get_range <- function(y_vals) {
    list(y_min = min(y_vals), y_max = max(y_vals))
}

plot_prepare <- function(name, g, x0, epsilon = 1e-6, max_iter = 100) {
    fixed_point <- x0

    # Intervallum pontosítás a fixpont algoritmus prototípus futtatásával
    for (i in 1:max_iter) {
        fixed_point <- g(fixed_point)
        if (abs(fixed_point - g(fixed_point)) < epsilon) break
    }
    
    zoom_factor <- abs(g(fixed_point) - fixed_point)
    if (zoom_factor < epsilon) {
        zoom_factor <- 1 
    }
    
    a <- fixed_point - zoom_factor
    b <- fixed_point + zoom_factor
    
    x_vals <- seq(a, b, length.out = 100)
    y_vals <- sapply(x_vals, g)
    range <- get_range(y_vals)
    
    plot(x_vals, x_vals, type = "l", col = "red", lwd = 2, xlab = "x", ylab = "y", 
         ylim = c(range$y_min, range$y_max), main = paste("Fixpont iteráció |", name))
    lines(x_vals, y_vals, col = "blue", lwd = 2)
    
    rect(a, range$y_min, b, range$y_max, border = "lightgray")
    
    return(range$y_min)
}

plot_iter <- function(iter, x, y, y_min) {
	lines(c(x, y), c(y, y), col = "black", lwd = 1)

	lines(c(x, x), c(x, y), col = "black", lwd = 1)

    if (iter == 1) {
        lines(c(x, x), c(y, y_min), col = "black", lwd = 1, lty = 2)
    }
    if (iter < 4) {
        lines(c(y, y), c(y, y_min), col = "black", lwd = 1, lty = 2)
    }
    if (iter < 5) {
        text(x, y_min, labels = paste0("x", iter - 1), pos = 1, col = "red", cex = 0.8)
    }
}

plot_finalize <- function(y, y_min) {
	points(y, y, col = "darkblue", pch = 18, lwd = 3, cex = 2.5)
	lines(c(y, y), c(y, y_min), col = "darkblue", lty = 2, lwd = 2)
    text(y, y_min, labels = "xm", pos = 1, col = "darkblue", cex = 0.8)

	legend("topleft", legend = c("y = x", "y = g(x)", "Iterációk", "Fixpont"), 
			col = c("red", "blue", "black", "darkblue"), 
			lty = c(1, 1, 1, 2), lwd = c(2, 2, 1, 2), pch = c(NA, NA, NA, 18), 
			pt.cex = c(NA, NA, NA, 2.5))
}

print_step <- function(iter, x, y) {
    cat("\nIteráció", iter, "\n")
    cat("x:", x, "\n")
    cat("y:", y, "\n")
}

# Fixpont iteráció algoritmus
fixpoint_iteration <- function(name, g, x0, epsilon = 1e-6, max_iter = 100) {
	y_min <- plot_prepare(name, g, x0, epsilon, max_iter)
    
	x <- x0
	x_prev <- x0
  
  	for (iter in 1:max_iter) {
		y <- g(x)
		
		plot_iter(iter, x, y, y_min)
		print_step(iter, x, y)

		if (abs(y - x) <= epsilon && abs(y - g(y)) <= epsilon) {
			plot_finalize(y, y_min)
			
			return(list(fixed_point = y, iterations = iter))
		}
		
		x_prev <- x
		x <- y
  	}
  
	cat("\nMaximum iteráció elérve. Nincs konvergencia.\n")
	return(list(fixed_point = x, iterations = max_iter))
}

read_choice <- function() {
    message("1: Lineáris konvergencia");
    message("2: Négyzetes konvergencia");
    message("3: Oszcillációs konvergencia");
    
    choice <- 1
    
    repeat {
        choice <- as.numeric(readline())
        
        if (!is.na(choice) && choice %in% 1:3) {
            return(choice)
        }
        
        cat("Érvénytelen érték, próbáld újra!\n")
    }
}

get_demo <- function(choice) {
    demos <- list(
        "1" = list( name = "Lineáris konvergencia",
                    g = function(x){ (x + 2) / 3 }, x0 = 0),
        "2" = list(name = "Négyzetes konvergencia",
                    g = function(x){ sqrt(x + 1) }, x0 = 1),
        "3" = list(name = "Oszcillációs konvergencia",
                    g = function(x){ cos(x) },      x0 = 1)
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
        cat("g(x) =", get_function_body(demo$g), "\n")
		cat("x0 =\t", demo$x0, "\n")

        result <- fixpoint_iteration(demo$name, demo$g, demo$x0)

        cat("\nKonvergencia elérve", result$iterations, "iteráció után.\n")
        cat("\nFixpont: ", result$fixed_point, "\n")

        message("\nÚjra futtatod? i / n")
        if (tolower(readline()) != "i") break
    }
}

main()
