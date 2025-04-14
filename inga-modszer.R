solve_block_tridiag <- function(A, B, C, F, verbose = FALSE) {
    M <- length(C)
    ALPHA <- vector("list", M)
    BETA <- vector("list", M)
    
    if (verbose) cat("\n=== ELŐREHALADÓ FÁZIS ===\n")
    ALPHA[[1]] <- C[[1]]
    BETA[[1]] <- F[[1]]
    
    for (i in 2:M) {
        if (verbose) cat("\nBlokk", i, "feldolgozása...\n")
        
        A_i     <- if (!is.null(A[[i]])) A[[i]]     else matrix(0, nrow=nrow(C[[1]]), ncol=ncol(C[[1]]))
        B_i_1   <- if (!is.null(B[[i-1]])) B[[i-1]] else matrix(0, nrow=nrow(C[[1]]), ncol=ncol(C[[1]]))
        
        inv_alpha_prev <- solve(ALPHA[[i-1]])
        if (verbose) {
            cat("α [",i-1,"] inverze:\n")
            print(round(inv_alpha_prev, 3))
        }
        
        ALPHA[[i]] <- C[[i]] - A_i %*% inv_alpha_prev %*% B_i_1
        BETA[[i]] <- F[[i]] - A_i %*% inv_alpha_prev %*% BETA[[i-1]]
        
        if (verbose) {
            cat("\nα [",i,"]:\n")
            print(round(ALPHA[[i]], 3))
            cat("\nβ [",i,"]:\n")
            print(round(BETA[[i]], 3))
        }
    }
    
    if (verbose) cat("\n=== VISSZAHALADÓ FÁZIS ===\n")
    X <- vector("list", M)
    X[[M]] <- solve(ALPHA[[M]]) %*% BETA[[M]]
    
    for (i in (M-1):1) {
        if (verbose) cat("\nBlokk", i, "feldolgozása...\n")
        B_i <- if (!is.null(B[[i]])) B[[i]] else matrix(0, nrow=nrow(C[[1]]), ncol=ncol(C[[1]]))
        
        X[[i]] <- solve(ALPHA[[i]]) %*% (BETA[[i]] - B_i %*% X[[i+1]])
        
        if (verbose) {
            cat("X [",i+1,"] értéke:\n")
            print(round(X[[i+1]], 3))
            cat("\nX [",i,"] végeredménye:\n")
            print(round(X[[i]], 3))
        }
    }
    
    return(X)
}

reconstruct_H <- function(A, B, C) {
    M <- length(C)
    N <- nrow(C[[1]])
    
    Hfull <- matrix(0, nrow=M*N, ncol=M*N)
    for (i in 1:M) {
        rowRange <- ((i-1)*N+1):(i*N)
        colRange <- rowRange
        Hfull[rowRange, colRange] <- C[[i]]
        
        if (i > 1) {
            rowRangeA <- ((i-1)*N+1):(i*N)
            colRangeA <- ((i-2)*N+1):((i-1)*N)
            Hfull[rowRangeA, colRangeA] <- A[[i]]
        }
        
        if (i < M) {
            rowRangeB <- ((i-1)*N+1):(i*N)
            colRangeB <- (i*N+1):((i+1)*N)
            Hfull[rowRangeB, colRangeB] <- B[[i]]
        }
    }
    return(Hfull)
}

read_choice <- function() {
    cat("\n=== INGA MÓDSZER PROGRAM === \n")
    message("1: Random mátrix")
    message("2: Demo mátrix")
    
    repeat {
        choice <- as.numeric(readline())
        if (!is.na(choice) && choice %in% 1:3) {
            return(choice)
        }
        cat("Érvénytelen választás, próbáld újra!\n")
    }
}

print_block_matrix <- function(A, B, C, M, N, demo_mode = FALSE) {
  cat("\n=== BLOKK TRIDIAGONÁLIS MÁTRIX STRUKTÚRA ===\n")
  for (i in 1:M) {
		line <- ""
		
		if (i > 1) {
			if (demo_mode) {
				line <- paste0(line, "A_", i, ":\n")
				print(A[[i]])
				cat("\n")
			} else {
				line <- paste0(line, "A_", i, " ")
			}
		}
		
		if (demo_mode) {
			line <- paste0(line, "C_", i, ":\n")
			print(C[[i]])
			cat("\n")
			} else {
			line <- paste0(line, "C_", i, " ")
		}
		
		if (i < M) {
			if (demo_mode) {
				line <- paste0(line, "B_", i, ":\n")
				print(B[[i]])
				cat("\n")
			} else {
				line <- paste0(line, "B_", i, " ")
			}
		}
		
		if (!demo_mode) {
			prefix <- if (i > 1) paste(rep("0 ", (i-2)*3), collapse = "") else ""
			suffix <- if (i < M) paste(rep("0 ", (M-i-1)*3), collapse = "")  else ""
			cat(prefix, line, suffix, "\n")
		}
	}
}

print_vectors <- function(F, M, demo_mode = FALSE) {
	cat("\n=== JOBB OLDALI VEKTOROK ===\n")
	for(i in 1:M) {
		if (demo_mode) {
		cat("F_", i, ":\n", round(F[[i]], 2), "\n\n")
		} else {
		cat("F_", i, "\n")
		}
	}
}

print_Hfull <- function(Hfull, M, N) {
	cat("\n=== RECONSTRUKÁLT TELJES MÁTRIX ===\n")
	if (M*N <= 6) {
		print(round(Hfull, 2))
	}
	else {
		cat("[", N, "x", N, " blokkokból épül fel, összes méret: ", M*N, "x", M*N, "]\n")
		cat("Példa blokk elrendezés:\n")
		print(Hfull[1:(2*N), 1:(3*N)])
	}
}

main <- function() {
    repeat {
        choice <- read_choice()

        if (choice == 1) {
            set.seed(123)
            M <- as.numeric(readline("Blokkok száma (M): "))
            N <- as.numeric(readline("Almátrix mérete (N): "))
            A <- vector("list", M)
            B <- vector("list", M)
            C <- vector("list", M)
            F <- vector("list", M)
            
            for (i in 1:M) {
                A[[i]] <- matrix(rnorm(N*N), nrow=N)
                B[[i]] <- matrix(rnorm(N*N), nrow=N)
                C[[i]] <- matrix(rnorm(N*N), nrow=N)
                F[[i]] <- rnorm(N)
            }
            
            A[[1]] <- matrix(0, nrow=N, ncol=N)
            B[[3]] <- matrix(0, nrow=N, ncol=N)

            print_block_matrix(A, B, C, M, N)
            print_vectors(F, M)
        }
        else if (choice == 2) {
            A <- list(
                matrix(0, 2, 2),
                matrix(c(1,0,0,1), 2),
                matrix(c(1,2,3,4), 2)
            )
            B <- list(
                matrix(c(2,1,1,2), 2),
                matrix(c(0,1,1,0), 2),
                matrix(0, 2, 2)
            )
            C <- list(
                matrix(c(4,1,1,4), 2),
                matrix(c(5,0,0,5), 2),
                matrix(c(6,1,1,6), 2)
            )
            F <- list(
                c(1, 2),
                c(3, 4),
                c(5, 6)
            )

            print_block_matrix(A, B, C, M=3, N=2, demo_mode=TRUE)
            print_vectors(F, M=3, demo_mode=TRUE)
        }
        
        X_solution <- solve_block_tridiag(A, B, C, F, verbose = TRUE)
        Hfull <- reconstruct_H(A, B, C)
        
        print_Hfull(Hfull, M, N)

        exact_solution <- solve(Hfull) %*% unlist(F)
        error <- norm(unlist(X_solution) - exact_solution)
        cat("Relatív hiba:", error/norm(exact_solution), "\n")
        
        message("\nÚjra futtatod? i / n")
        if (tolower(readline()) != "i") break
    }
}

main()