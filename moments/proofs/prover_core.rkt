#lang racket

(require racket/trace)
(current-prefix-in " ")
(current-prefix-out " ")

(provide symbolic-diff
         symbolic-simp-rule
         symbolic-simp
         is-real
         flux-deriv-replace
         symbolic-roe-function
         is-non-negative
         variable-transform
         symbolic-simp-positive-rule
         symbolic-simp-positive
         evaluate-limit-rule
         evaluate-limit
         prove-lax-friedrichs-scalar-1d-hyperbolicity
         prove-lax-friedrichs-scalar-1d-cfl-stability
         prove-lax-friedrichs-scalar-1d-local-lipschitz
         prove-roe-scalar-1d-hyperbolicity
         prove-roe-scalar-1d-flux-conservation
         prove-flux-limiter-symmetry
         prove-flux-limiter-tvd)

;; Lightweight symbolic differentiator (differentiates expr with respect to var).
(define (symbolic-diff expr var)
  (match expr
    ;; If expr is a symbol, then it either differentiates to 1 (if it's equal to var), or 0 otherwise.
    [(? symbol? symb) (cond
                        [(eq? symb var) 1.0]
                        [else 0.0])]

    ;; If expr is a numerical constant, then it differentiates to 0.
    [(? number?) 0.0]

    ;; If expr is a sum of the form (+ expr1 expr2 ...), then it differentiates to a sum of derivatives (+ expr1' expr2' ...), by linearity.
    [`(+ . ,terms)
     `(+ ,@(map (lambda (term) (symbolic-diff term var)) terms))]
    ;; Likewise for differences of the form (- expr1 expr2 ...), which differentiate to (- expr1' expr2' ...), by linearity.
    [`(- . ,terms)
     `(- ,@(map (lambda (term) (symbolic-diff term var)) terms))]

    ;; If expr is a product of the form (* expr1 expr2 ...), then it differentiates to (+ (* expr1' expr2 ...) (* expr1 expr2' ...) ...), by the product rule.
    [`(* . ,terms)
     (define n (length terms))
     (define (mult xs) (cons '* xs)) ; Multiplication helper function.

     ((lambda (sums) (cond
                      [(null? (cdr sums)) (car sums)]
                      [else (cons '+ sums)]))
      (let loop ([i 0])
        (cond
          [(= i n) `()]
          [else
           ;; Evaluate the derivative of the i-th term in the product.
           (let ([di (symbolic-diff (list-ref terms i) var)])
             (cons
              (mult (for/list ([j (in-range n)])
                      (cond
                        [(= j i) di]
                        [else (list-ref terms j)])))
              (loop (add1 i))))])))]

    ;; If expr is a quotient of the form (/ expr1 expr2), then it differentiates to (/ (- (* expr2 expr1') (expr1 expr2') (* expr2 expr2)), by the quotient rule.
    [`(/ ,x ,y)
     `(/ (- (* ,y ,(symbolic-diff x var)) (* ,x ,(symbolic-diff y var))) (* ,y ,y))]

    ;; If expr is an absolute value of the form (abs expr1), then it differentiates to (sgn expr1').
    [`(abs ,arg)
     `(* (sgn ,arg) ,(symbolic-diff arg var))]

    ;; If expr is a sign function of the form (sgn expr1), then it differentiates to 0.0.
    [`(sgn ,arg) 0.0]
    
    ;; Otherwise, return false.
    [else #f]))

;; Lightweight symbolic simplification rules (simplifies expr using only correctness-preserving algebraic transformations).
(define (symbolic-simp-rule expr)
  (match expr    
    ;; If expr is of the form (0 + x) or (0.0 + x), then simplify to x.
    [`(+ 0 ,x) `,x]
    [`(+ 0.0 ,x) `,x]
    [`(+ -0.0 ,x) `,x]

    ;; If expr is of the form (1 * x) or (1.0 * x), then simplify to x.
    [`(* 1 ,x) `,x]
    [`(* 1.0 ,x) `,x]

    ;; If expr is of the form (0 * x) or (0.0 * x), then simplify to 0 or 0.0.
    [`(* 0 ,x) 0]
    [`(* 0.0 ,x) 0.0]
    [`(* -0.0 ,x) 0.0]

    ;; If expr is of the form (x - 0) or (x - 0.0), then simplify to x.
    [`(- ,x 0) `,x]
    [`(- ,x 0.0) `,x]
    [`(- ,x -0.0) `,x]

    ;; If expr is of the form (0 - x) or (0.0 - x), then simplify to (-1 * x) or (-1.0 * x).
    [`(- 0 ,x) `(* -1 ,x)]
    [`(- 0.0 ,x) `(* -1.0 ,x)]
    [`(- -0.0 ,x) `(* -1.0 ,x)]

    ;; If expr is of the form (x / 1) or (x / 1.0), then simplify to x.
    [`(/ ,x 1) `,x]
    [`(/ ,x 1.0) `,x]

    ;; Enforce right associativity of addition: if expr is of the form ((x + y) + z) or (x + y + z), then simplify to (x + (y + z)).
    [`(+ (+ ,x ,y) ,z) `(+ ,x (+ ,y ,z))]
    [`(+ ,x ,y ,z) `(+ (+ ,x ,y) ,z)]

    ;; Enforce right associativity of multiplication: if expr is of the form ((x * y) * z) or (x * y * z), then simplify to (x * (y * z)).
    [`(* (* ,x ,y) ,z) `(* ,x (* ,y ,z))]
    [`(* ,x ,y ,z) `(* (* ,x ,y) ,z)]

    ;; If expr is of the form (x + y) for numeric x and y, then just evaluate the sum. Likewise for differences.
    [`(+ ,(and x (? number?)) ,(and y (? number?))) (+ x y)]
    [`(- ,(and x (? number?)) ,(and y (? number?))) (- x y)]

    ;; If expr is of the form (x * y) for numeric x and y, then just evaluate the product. Likewise for quotients
    [`(* ,(and x (? number?)) ,(and y (? number?))) (* x y)]
    [`(/ ,(and x (? number?)) ,(and y (? number?))) (/ x y)]

    ;; If expr is of the form (x * (y + z)) for numeric x, y and z, then just evaluate the product and sum.
    [`(* ,(and x (? number?)) (+ ,(and y (? number?)) ,(and z (? number?)))) (* x (+ y z))]

    ;; Enforce (reverse) distributive property: if expr is of the form ((a * x) + (b * x)), then simplify to ((a + b) * x).
    [`(+ (* ,a, x) (* ,b ,x)) `(* (+ ,a ,b) ,x)]

    ;; If expr is of the form (x * (y * z)) for numeric numeric x and y, then evaluate the product of x and y.
    [`(* ,(and x (? number?)) (* ,(and y (? number?)) ,z)) `(* ,(* x y) ,z)]

    ;; Move numbers to the left: if expr is of the form (x + y) for non-numeric x but numeric y, then simplify to (y + x).
    [`(+ ,(and x (not (? number?))) ,(and y (? number?))) `(+ ,y ,x)]

    ;; Move numbers to the left: if expr is of the form (x * y) for non-numeric x but numeric y, then simplify to (y * x).
    [`(* ,(and x (not (? number?))) ,(and y (? number?))) `(* ,y ,x)]

    ;; If expr is of the form sqrt(x * x) or (sqrt(x) * sqrt(x)), then simplify to x.
    [`(sqrt (* ,x ,x)) `,x]
    [`(* (sqrt ,x) (sqrt ,x)) `,x]

    ;; If expr is of the form (sqrt(x) * (y * sqrt(x))), then simplify to (y * x).
    [`(* (sqrt,x) (* ,y (sqrt ,x))) `(* ,y ,x)]
    ;; Likewise, if expr is of the form (sqrt(x) * (sqrt(x) * y)), then simplify to (x * y).
    [`(* (sqrt,x) (* (sqrt ,x) ,y)) `(* ,x ,y)]

    ;; If expr is of the form sqrt(x * y), then simplify to (sqrt(x) * sqrt(y)).
    [`(sqrt (* ,x ,y)) `(* (sqrt ,x) (sqrt ,y))]

    ;; If expr if of the form sqrt(x) for numeric x, then just evaluate the square root.
    [`(sqrt ,(and x (? number?))) (sqrt x)]

    ;; If expr is of the form max(x, y) or min(x, y) for numeric x and y, then just evaluate the maximum/minimum.
    [`(max ,(and x (? number?)) ,(and y (? number?))) (max x y)]
    [`(min ,(and x (? number?)) ,(and y (? number?))) (min x y)]

    ;; If expr is of the form abs(x) for numeric x, then just evaluate the absolute value.,
    [`(abs ,(and x (? number?))) (abs x)]

    ;; If expr is of the form abs(-1 * x) or abs(-1.0 * x), then simplify to abs(x).
    [`(abs (* -1 ,x)) `(abs ,x)]
    [`(abs (* -1.0 ,x)) `(abs ,x)]

    ;; If expr is of the form (0 - (x * y)) or (0.0 - (x * y)), then simplify to ((0 - x) * y) or ((0.0 - x) * y).
    [`(- 0 (* ,x ,y)) `(* (- 0 ,x) ,y)]
    [`(- 0.0 (* ,x ,y)) `(* (- 0.0 ,x) ,y)]
    [`(- -0.0 (* ,x ,y)) `(* (- 0.0 ,x) ,y)]

    ;; If expr is of the form (x + x), thens implify to (2.0 * x).
    [`(+ ,x ,x) `(* 2.0 ,x)]

    ;; If expr is of the form ((x * y) / (x * z)), then simplify to (y / z).
    [`(/ (* ,x ,y) (* ,x ,z)) `(/ ,y ,z)]

    ;; If expr is of the form ((x / y) * (x / y)), then simplify to ((x * x) / (y * y)).
    [`(* (/ ,x ,y) (/ ,x ,y)) `(/ (* ,x ,x) (* ,y ,y))]

    ;; If expr is of the form (x * (y * z)) for numeric y and non-numeric x and z, then simplify to (y * (x * z)).
    [`(* ,(and x (not (? number?))) (* ,(and y (? number?)) ,(and z (not (? number?))))) `(* ,y (* ,x ,z))]

    ;; Enforce distributive property: if expr is of the form (x * (a + b)), then simplify to ((x * a) + (x * b)).
    [`(* ,x (+ ,a ,b)) `(+ (* ,x ,a) (* ,x ,b))]

    ;; If expr is of the form (x * (-y / z)), then simplify to (-x * (y / z)).
    [`(* ,x (/ (* -1 ,y) ,z)) `(* (* -1 ,x) (/ ,y ,z))]
    [`(* ,x (/ (* -1.0 ,y) ,z)) `(* (* -1.0 ,x ) (/ ,y ,z))]

    ;; If expr is of the form ((x * y) / z) for numeric x, then simplify to (x * (y / z)).
    [`(/ (* ,(and x (? number?)) ,y) ,z) `(* ,x (/ ,y ,z))]

    ;; If expr is of the form ((a * x) + (y + (b * x))) for numeric a and b, then simplify to (((a + b) * x) + y).
    [`(+ (* ,(and a (? number?)) ,x) (+ ,y (* ,(and b (? number?)) ,x))) `(+ (* (+ ,a ,b) ,x) ,y)]

    ;; If expr is of the form (a + (x / y)) or (-a + (x / y)) for symbolic a, then simplify to ((x / y) + a) or ((x / y) - a).
    [`(+ ,(and a (? symbol?)) (/ ,x ,y)) `(+ (/ ,x ,y) ,a)]
    [`(+ (* -1 ,(and a (? symbol?))) (/ ,x ,y)) `(- (/ ,x ,y) ,a)]
    [`(+ (* -1.0 ,(and a (? symbol?))) (/ ,x ,y)) `(- (/ ,x ,y) ,a)]

    ;; Enforce (reverse) distributive property: if expr is of the form ((a * x) - (a * y)), then simplify to (a * (x - y)).
    [`(- (* ,a ,x) (* ,a ,y)) `(* ,a (- ,x ,y))]

    ;; If expr is of the form (((a * x) + (a * y)) * (x - y)), then simplify to ((a * (x * x)) - (a * (y * y))).
    [`(* (+ (* ,a ,x) (* ,a ,y)) (- ,x ,y)) `(- (* ,a (* ,x ,x)) (* ,a (* ,y ,y)))]

    ;; If expr is of the form (0 / x) or (0.0 / x), then simplify to 0 or 0.0.
    [`(/ 0 ,x) 0]
    [`(/ 0.0 ,x) 0.0]
    [`(/ -0.0 ,x) 0.0]

    ;; If expr is of the form (x / x), then simplify to 1.0
    [`(/ ,x ,x) 1.0]

    ;; If expr is of the form (x * (y / z)) for numeric x and y, then evaluate the product to yield ((x * y) / z).
    [`(* ,(and x (? number?)) (/ ,(and y (? number?)) ,z)) `(/ ,(* x y) ,z)]
    ;; Likewise, if expr is of the form ((x / y) / z) for numeric x and z, then evaluate the quotient to yield ((x / z) / y).
    [`(/ (/ ,(and x (? number?)) ,y) ,(and z (? number?))) `(/ ,(/ x z) ,y)]

    ;; If expr is of the form ((x / y) / x), then simplify to (1.0 / y).
    [`(/ (/ ,x ,y) ,x) `(/ 1.0 ,y)]

    ;; If expr is of the form ((x / y) / (z + (x / y))), or ((x / y) / ((x / y) + z), then simplify to (x / ((z * y) + x)) or (x / (x + (z * y))).
    [`(/ (/ ,x ,y) (+ ,z (/ ,x ,y))) `(/ ,x (+ (* ,z ,y) ,x))]
    [`(/ (/ ,x ,y) (+ (/ ,x ,y) ,z)) `(/ ,x (+ ,x (* ,z ,y)))]

    ;; If expr is of the form ((x + y) / z) or ((x - y) / z), then simplify to ((x / z) + (y / z)) or ((x / z) - (y / z)).
    [`(/ (+ ,x ,y) ,z) `(+ (/ ,x ,z) (/ ,y ,z))]
    [`(/ (- ,x ,y) ,z) `(- (/ ,x ,z) (/ ,y ,z))]

    ;; If expr is a sum of the form (x + y + ...), then apply symbolic simplification to each term x, y, ... in the sum.
    [`(+ . ,terms)
     `(+ ,@(map (lambda (term) (symbolic-simp-rule term)) terms))]
    ;; Likewise for differences.
    [`(- . ,terms)
     `(- ,@(map (lambda (term) (symbolic-simp-rule term)) terms))]

    ;; If expr is a product of the form (x * y * ...), then apply symbolic simplification to each term x, y, ... in the product.
    [`(* . ,terms)
     `(* ,@(map (lambda (term) (symbolic-simp-rule term)) terms))]
    ;; Likewise for quotients.
    [`(/ . ,terms)
     `(/ ,@(map (lambda (term) (symbolic-simp-rule term)) terms))]

    ;; If expr is of the form sqrt(expr1), then apply symbolic simplification to the interior expr1.
    [`(sqrt ,arg)
     `(sqrt ,(symbolic-simp-rule arg))]

    ;; If expr is of the form abs(expr1), then apply symbolic simplification to the interior expr1.
    [`(abs ,arg)
     `(abs ,(symbolic-simp-rule arg))]

    ;; If expr is of the form max(x, y, z) or min(x, y, z), then simplify to max(max(x, y), z) or min(min(x, y), z).
    [`(max ,x ,y ,z) `(max (max ,x ,y) ,z)]
    [`(min ,x ,y ,z) `(min (min ,x ,y) ,z)]

    ;; If expr is of the form max(x, y), then simplify to ((0.5 * (x + y)) + (0.5 * abs(x - y))).
    [`(max ,x ,y) `(+ (* 0.5 (+ ,x ,y)) (* 0.5 (abs (- ,x ,y))))]

    ;; If expr is of the form min(x, y), then simplify to ((0.5 * (x + y)) - (0.5 * abs(x - y))).
    [`(min ,x ,y) `(- (* 0.5 (+ ,x ,y)) (* 0.5 (abs (- ,x ,y))))]

    ;; If expr is a complex number whose imaginary part is equal to 0.0 or -0.0, then simplify to Re(expr).
    [(? (lambda (arg)
         (and (number? arg) (not (real? arg )) (equal? (imag-part arg) 0.0)))) (real-part expr)]
    [(? (lambda (arg)
         (and (number? arg) (not (real? arg )) (equal? (imag-part arg) -0.0)))) (real-part expr)]

    ;; Otherwise, return the expression.
    [else expr]))

;; Recursively apply the symbolic simplification rules until the expression stops changing (fixed point).
(define (symbolic-simp expr)
  (cond
    [(equal? (symbolic-simp-rule expr) expr) expr]
    [else (symbolic-simp (symbolic-simp-rule expr))]))

;; Recursively determine whether an expression corresponds to a real number.
(define (is-real expr cons-vars parameters)
  (match expr
    ;; Real numbers are trivially real.
    [(? real?) #t]

    ;; Conserved variables are assumed to be real (this is enforced elsewhere).
    [(? (lambda (arg)
          (not (equal? (member arg cons-vars) #f)))) #t]

    ;; Simulation parameters are assumed to be real (this is enforced elsewhere).
    [(? (lambda (arg)
          (and (not (empty? parameters)) (ormap (lambda (parameter)
                                                  (equal? arg (list-ref parameter 1))) parameters)))) #t]

    ;; The outcome of a conditional operation is real if both branches yield real numbers.
    [`(cond
        [,cond1 ,expr1]
        [else ,expr2])
     (and (is-real expr1 cons-vars parameters) (is-real expr2 cons-vars parameters))]

    ;; The sum, difference, product, or quotient of two real numbers is always real.
    [`(+ . ,terms)
     (andmap (lambda (term) (is-real term cons-vars parameters)) terms)]
    [`(- . ,terms)
     (andmap (lambda (term) (is-real term cons-vars parameters)) terms)]
    [`(* . ,terms)
     (andmap (lambda (term) (is-real term cons-vars parameters)) terms)]
    [`(/ . ,terms)
     (andmap (lambda (term) (is-real term cons-vars parameters)) terms)]

    ;; Otherwise, assume false.
    [else #f]))

;; Recursively replace conserved variable expressions within the flux derivative expression (for Roe functions).
(define (flux-deriv-replace flux-deriv-expr cons-expr new-cons-expr)
  (match flux-deriv-expr
    ;; If the flux derivative expression is just the conserved variable expression, then return the new conserved variable expression.
    [(? (lambda (arg)
          (equal? arg cons-expr))) new-cons-expr]

    ;; If the flux derivative expression consists of a sum, difference, product, or quotient, then recursively apply replacement to each term.
    [`(+ . ,terms)
     `(+ ,@(map (lambda (term) (flux-deriv-replace term cons-expr new-cons-expr)) terms))]
    [`(- . ,terms)
     `(- ,@(map (lambda (term) (flux-deriv-replace term cons-expr new-cons-expr)) terms))]
    [`(* . ,terms)
     `(* ,@(map (lambda (term) (flux-deriv-replace term cons-expr new-cons-expr)) terms))]
    [`(/ . ,terms)
     `(/ ,@(map (lambda (term) (flux-deriv-replace term cons-expr new-cons-expr)) terms))]

    ;; Otherwise, return the flux derivative expression.
    [else flux-deriv-expr]))

;; Compute the symbolic Roe function (averaged flux derivative).
(define (symbolic-roe-function flux-deriv-expr cons-expr)
  (symbolic-simp `(+ (* 0.5 ,(flux-deriv-replace flux-deriv-expr cons-expr (string->symbol (string-append (symbol->string cons-expr) "L"))))
                     (* 0.5 ,(flux-deriv-replace flux-deriv-expr cons-expr (string->symbol (string-append (symbol->string cons-expr) "R")))))))

;; Determine whether an expression is non-negative.
(define (is-non-negative expr parameters)
  (match expr
    ;; A non-negative number is, trivially, non-negative.
    [(? (lambda (arg)
          (and (number? arg) (or (>= arg 0) (>= arg 0.0))))) #t]

    ;; Simulation parameters that are non-negative are, trivially, non-negative.
    [(? (lambda (arg)
          (and (not (empty? parameters)) (ormap (lambda (parameter)
                                                  (and (equal? arg (list-ref parameter 1))
                                                       (or (>= (list-ref parameter 2) 0)
                                                           (>= (list-ref parameter 2) 0.0)))) parameters)))) #t]

    ;; The sum, product, or quotient of two non-negative numbers is always non-negative.
    [`(+ ,x ,y) (and (is-non-negative x parameters) (is-non-negative y parameters))]
    [`(* ,x ,y) (and (is-non-negative x parameters) (is-non-negative y parameters))]
    [`(/ ,x ,y) (and (is-non-negative x parameters) (is-non-negative y parameters))]

    ;; Otherwise, assume false.
    [else #f]))

;; Recursively ransform all occurrences of a given variable within an expression to a new variable.
(define (variable-transform expr var new-var)
  (cond
    ;; Replace any occurrence of var in expr with new-var.
    [(symbol? expr) (cond
                      [(equal? expr var) new-var]
                      [else expr])]

    ;; Recursively apply variable-transform to all subexpressions.
    [(pair? expr) (map (lambda (subexpr)
                         (variable-transform subexpr var new-var)) expr)]

    ;; Otherwise, return the expression.
    [else expr]))

;; Lightweight symbolic simplification rules, assuming strict positivity of pos-var.
(define (symbolic-simp-positive-rule expr pos-var)
  (match expr
    ;; If expr is of the form (abs(x) / y), with y strictly positive, then simplify to abs(x / y).
    [`(/ (abs ,x) ,pos-var) `(abs (/ ,x ,pos-var))]

    ;; If expr is of the form abs(x), abs(1 / x) or abs(1.0 / x), with x strictly positive, then simplify to x, (1 / x) or (1.0 / x).
    [`(abs ,pos-var) pos-var]
    [`(abs (/ 1 ,pos-var)) `(/ 1 ,pos-var)]
    [`(abs (/ 1.0 ,pos-var)) `(/ 1.0 ,pos-var)]

    ;; If expr is of the form sgn(x), with x strictly positive, then simplify to 1.0.
    [`(sgn ,pos-var) 1.0]

    ;; If expr is of the form abs(expr1), then apply symbolic simplification to the interior expr1.
    [`(abs ,x)
     `(abs ,(symbolic-simp-positive-rule x pos-var))]

    ;; If expr is of the form (max(x, y) / z) or (min(x, y) / z), with z strictly positive, then simplify to max((x / z), (y / z)) or min((x / z), (y / z)).
    [`(/ (max ,x ,y) ,pos-var) `(max (/ ,y ,pos-var) (/ ,x ,pos-var))]
    [`(/ (min ,x ,y) ,pos-var) `(min (/ ,y ,pos-var) (/ ,x ,pos-var))]

    ;; If expr is of the form (max(x, y, z) / w) or (min(x, y, z) / w), with w strictly positive, then simplify to max((x / w), (y / w), (z / w)) or min((x / w), (y / w), (z / w)).
    [`(/ (max ,x ,y ,z) ,pos-var) `(max (/ ,z ,pos-var) (/ ,y ,pos-var) (/ ,x ,pos-var))]
    [`(/ (min ,x ,y, z) ,pos-var) `(min (/ ,z ,pos-var) (/ ,y ,pos-var) (/ ,x ,pos-var))]

    ;; If expr is a max or a min of the form max(x, y, ...) or min(x, y, ...), then apply symbolic simplification to each term x, y, ... in the function.
    [`(max . ,terms)
     `(max ,@(map (lambda (term) (symbolic-simp-positive-rule term pos-var)) terms))]
    [`(min . ,terms)
     `(min ,@(map (lambda (term) (symbolic-simp-positive-rule term pos-var)) terms))]

    ;; If expr is a sum of the form (x + y + ...), then apply symbolic simplification to each term x, y, ... in the sum.
    [`(+ . ,terms)
     `(+ ,@(map (lambda (term) (symbolic-simp-positive-rule term pos-var)) terms))]
    ;; Likewise for differences.
    [`(- . ,terms)
     `(- ,@(map (lambda (term) (symbolic-simp-positive-rule term pos-var)) terms))]

    ;; If expr is a product of the form (x * y * ...), then apply symbolic simplification to each term x, y, ... in the product.
    [`(* . ,terms)
     `(* ,@(map (lambda (term) (symbolic-simp-positive-rule term pos-var)) terms))]
    ;; Likewise for quotients.
    [`(/ . ,terms)
     `(/ ,@(map (lambda (term) (symbolic-simp-positive-rule term pos-var)) terms))]

    ;; Otherwise, return the expression.
    [else expr]))

;; Recursively apply the symbolic simplification rules (assuming strict positivity of pos-var) until the expression stops changing (fixed point).
(define (symbolic-simp-positive expr pos-var)
  (cond
    [(equal? (symbolic-simp-positive-rule expr pos-var) expr) expr]
    [else (symbolic-simp-positive (symbolic-simp-positive-rule expr pos-var) pos-var)]))

;; Lightweight symbolic limit evaluation rules (computes limit of expr as var approaches lim).
(define (evaluate-limit-rule expr var lim)
  (match expr
    ;; If expr is of the form max(x, y) for numeric x and y, then just evaluate the maximum of the pair. Likewise for minima.
    [`(max ,(and x (? number?)) ,(and y (? number?))) (max x y)]
    [`(min ,(and x (? number?)) ,(and y (? number?))) (min x y)]

    ;; If expr is of the form max(x, y, z) for numeric x, y and z, then just evaluate the maximum of the triple. Likewise for minima.
    [`(max ,(and x (? number?)) ,(and y (? number?)) ,(and z (? number?))) (max x y z)]
    [`(min ,(and x (? number?)) ,(and y (? number?)) ,(and z (? number?))) (min x y z)]

    ;; If expr is of the form abs(x) for numeric x, then just evaluate the absolute value.
    [`(abs ,(and x (? number?))) (abs x)]

    ;; If expr is of the form (x + y) for numeric x and y, then just evaluate the sum. Likewise for differences.
    [`(+ ,(and x (? number?)) ,(and y (? number?))) (+ x y)]
    [`(- ,(and x (? number?)) ,(and y (? number?))) (- x y)]

    ;; If expr is of the form (x * y) for numeric x and y, then just evaluate the product. Likewise for quotients.
    [`(* ,(and x (? number?)) ,(and y (? number?))) (* x y)]
    [`(/ ,(and x (? number?)) ,(and y (? number?))) (/ x y)]

    ;; If expr is of the form max(expr1, expr2), then evaluate the limits of the interior expr1 and expr2. Likewise for minima.
    [`(max ,x ,y) `(max ,(evaluate-limit-rule x var lim) ,(evaluate-limit-rule y var lim))]
    [`(min ,x ,y) `(min ,(evaluate-limit-rule x var lim) ,(evaluate-limit-rule y var lim))]

    ;; If expr is of the form max(expr1, expr2, expr3), then evaluate the limits of the interior expr1, expr2 and expr3. Likewise for minima.
    [`(max ,x ,y ,z) `(max ,(evaluate-limit-rule x var lim) ,(evaluate-limit-rule y var lim) ,(evaluate-limit-rule z var lim))]
    [`(min ,x ,y ,z) `(min ,(evaluate-limit-rule x var lim) ,(evaluate-limit-rule y var lim) ,(evaluate-limit-rule z var lim))]

    ;; If expr is of the form abs(expr1), then evaluate the limit of the interior expr1.
    [`(abs ,x) `(abs ,(evaluate-limit-rule x var lim))]

    ;; If expr is a sum of the form (x + y + ...), then evaluate the limits each term x, y, ... in the sum.
    [`(+ . ,terms)
     `(+ ,@(map (lambda (term) (evaluate-limit-rule term var lim)) terms))]
    ;; Likewise for differences.
    [`(- . ,terms)
     `(- ,@(map (lambda (term) (evaluate-limit-rule term var lim)) terms))]

    ;; If expr is a product of the form (x * y * ...), then evaluate the limits each term x, y, ... in the product.
    [`(* . ,terms)
     `(* ,@(map (lambda (term) (evaluate-limit-rule term var lim)) terms))]
    ;; Likewise for quotients.
    [`(/ . ,terms)
     `(/ ,@(map (lambda (term) (evaluate-limit-rule term var lim)) terms))]

    ;; Otherwise, return the expression.
    [else expr]))

;; Recursively apply the limit evaluation rules until the expression stops changing (fixed point).
(define (evaluate-limit expr var limit)
  (define limit-val (variable-transform expr var limit))
  
  (cond
    [(equal? (evaluate-limit-rule limit-val var limit) expr) expr]
    [else (evaluate-limit (evaluate-limit-rule limit-val var limit) var limit)]))

;; ----------------------------------------------------------------------------------------
;; Prove hyperbolicity of the Lax–Friedrichs (Finite-Difference) Solver for a 1D Scalar PDE
;; ----------------------------------------------------------------------------------------
(define (prove-lax-friedrichs-scalar-1d-hyperbolicity pde
                                                      #:nx [nx 200]
                                                      #:x0 [x0 0.0]
                                                      #:x1 [x1 2.0]
                                                      #:t-final [t-final 1.0]
                                                      #:cfl [cfl 0.95]
                                                      #:init-func [init-func `(cond
                                                                                [(< x 1.0) 1.0]
                                                                                [else 0.0])])
   "Prove that the Lax-Friedrichs finite-difference method preserves hyperbolicity for the 1D scalar PDE specified by `pde`. 
  - `nx` : Number of spatial cells.
  - `x0`, `x1` : Domain boundaries.
  - `t-final`: Final time.
  - `cfl`: CFL coefficient.
  - `init-func`: Racket expression for the initial condition, e.g. piecewise constant."

  (define cons-expr (hash-ref pde 'cons-expr))
  (define flux-expr (hash-ref pde 'flux-expr))
  (define parameters (hash-ref pde 'parameters))

  (trace is-real)
  (trace symbolic-simp)
  (trace symbolic-simp-rule)
  (trace symbolic-diff)
  
  (define out (cond
    ;; Check whether the CFL coefficient is greater than 0 and less than or equal to 1 (otherwise, return false).
    [(or (<= cfl 0) (> cfl 1)) #f]
    
    ;; Check whether the number of spatial cells is at least 1 and the right domain boundary is set to the right of the left boundary (otherwise, return false)
    [(or (< nx 1) (>= x0 x1)) #f]
    
    ;; Check whether the final simulation time is non-negative (otherwise, return false).
    [(< t-final 0) #f]

    ;; Check whether the simulation parameter(s) correspond to real numbers (otherwise, return false).
    [(not (or (empty? parameters) (andmap (lambda (parameter)
                                            (is-real (list-ref parameter 2) (list cons-expr) parameters)) parameters))) #f]

    ;; Check whether the initial condition(s) correspond to real numbers (otherwise, return false).
    [(not (is-real init-func (list cons-expr) parameters)) #f]
    
    ;; Check whether the derivative of the flux function is real (otherwise, return false).
    [(not (is-real (symbolic-simp (symbolic-diff flux-expr cons-expr)) (list cons-expr) parameters)) #f]

    ;; Otherwise, return true.
    [else #t]))

  (untrace is-real)
  (untrace symbolic-simp)
  (untrace symbolic-simp-rule)
  (untrace symbolic-diff)
  
  out)
(trace prove-lax-friedrichs-scalar-1d-hyperbolicity)

;; ----------------------------------------------------------------------------------------
;; Prove CFL stability of the Lax–Friedrichs (Finite-Difference) Solver for a 1D Scalar PDE
;; ----------------------------------------------------------------------------------------
(define (prove-lax-friedrichs-scalar-1d-cfl-stability pde
                                                      #:nx [nx 200]
                                                      #:x0 [x0 0.0]
                                                      #:x1 [x1 2.0]
                                                      #:t-final [t-final 1.0]
                                                      #:cfl [cfl 0.95]
                                                      #:init-func [init-func `(cond
                                                                                [(< x 1.0) 1.0]
                                                                                [else 0.0])])
   "Prove that the Lax-Friedrichs finite-difference method is CFL stable for the 1D scalar PDE specified by `pde`. 
  - `nx` : Number of spatial cells.
  - `x0`, `x1` : Domain boundaries.
  - `t-final`: Final time.
  - `cfl`: CFL coefficient.
  - `init-func`: Racket expression for the initial condition, e.g. piecewise constant."

  (define cons-expr (hash-ref pde 'cons-expr))
  (define flux-expr (hash-ref pde 'flux-expr))
  (define max-speed-expr (hash-ref pde 'max-speed-expr))
  (define parameters (hash-ref pde 'parameters))

  (trace is-real)
  (trace symbolic-simp)
  (trace symbolic-simp-rule)
  (trace symbolic-diff)

  (define out (cond
    ;; Check whether the CFL coefficient is greater than 0 and less than or equal to 1 (otherwise, return false).
    [(or (<= cfl 0) (> cfl 1)) #f]
    
    ;; Check whether the number of spatial cells is at least 1 and the right domain boundary is set to the right of the left boundary (otherwise, return false)
    [(or (< nx 1) (>= x0 x1)) #f]
    
    ;; Check whether the final simulation time is non-negative (otherwise, return false).
    [(< t-final 0) #f]

    ;; Check whether the simulation parameter(s) correspond to real numbers (otherwise, return false).
    [(not (or (empty? parameters) (andmap (lambda (parameter)
                                            (is-real (list-ref parameter 2) (list cons-expr) parameters)) parameters))) #f]

    ;; Check whether the initial condition(s) correspond to real numbers (otherwise, return false).
    [(not (is-real init-func (list cons-expr) parameters)) #f]
    
    ;; Check whether the absolute value of the derivative of the flux function is symbolically equivalent to the maximum wave-speed estimate (otherwise, return false).
    [(not (equal? (symbolic-simp `(abs ,(symbolic-diff flux-expr cons-expr)))
                  (symbolic-simp max-speed-expr))) #f]

    ;; Otherwise, return true.
    [else #t]))

  (untrace is-real)
  (untrace symbolic-simp)
  (untrace symbolic-simp-rule)
  (untrace symbolic-diff)

  out)
(trace prove-lax-friedrichs-scalar-1d-cfl-stability)

;; ------------------------------------------------------------------------------------------------------------------------------------
;; Prove local Lipschitz continuity of the discrete flux function for the Lax–Friedrichs (Finite-Difference) Solver for a 1D Scalar PDE
;; ------------------------------------------------------------------------------------------------------------------------------------
(define (prove-lax-friedrichs-scalar-1d-local-lipschitz pde
                                                        #:nx [nx 200]
                                                        #:x0 [x0 0.0]
                                                        #:x1 [x1 2.0]
                                                        #:t-final [t-final 1.0]
                                                        #:cfl [cfl 0.95]
                                                        #:init-func [init-func `(cond
                                                                                  [(< x 1.0) 1.0]
                                                                                  [else 0.0])])
   "Prove that the Lax-Friedrichs finite-difference method has a discrete flux function that satisfies local Lipschitz continuity for the 1D scalar PDE specified by `pde`. 
  - `nx` : Number of spatial cells.
  - `x0`, `x1` : Domain boundaries.
  - `t-final`: Final time.
  - `cfl`: CFL coefficient.
  - `init-func`: Racket expression for the initial condition, e.g. piecewise constant."

  (define cons-expr (hash-ref pde 'cons-expr))
  (define flux-expr (hash-ref pde 'flux-expr))
  (define parameters (hash-ref pde 'parameters))

  (trace is-real)
  (trace symbolic-simp)
  (trace symbolic-simp-rule)
  (trace symbolic-diff)
  (trace is-non-negative)

  (define out (cond
    ;; Check whether the CFL coefficient is greater than 0 and less than or equal to 1 (otherwise, return false).
    [(or (<= cfl 0) (> cfl 1)) #f]
    
    ;; Check whether the number of spatial cells is at least 1 and the right domain boundary is set to the right of the left boundary (otherwise, return false)
    [(or (< nx 1) (>= x0 x1)) #f]
    
    ;; Check whether the final simulation time is non-negative (otherwise, return false).
    [(< t-final 0) #f]

    ;; Check whether the simulation parameter(s) correspond to real numbers (otherwise, return false).
    [(not (or (empty? parameters) (andmap (lambda (parameter)
                                            (is-real (list-ref parameter 2) (list cons-expr) parameters)) parameters))) #f]

    ;; Check whether the initial condition(s) correspond to real numbers (otherwise, return false).
    [(not (is-real init-func (list cons-expr) parameters)) #f]
    
    ;; Check whether the flux function is convex, i.e. that the second derivative of the flux function is strictly non-negative (otherwise, return false).
    [(let ([deriv (symbolic-simp (symbolic-diff (symbolic-simp (symbolic-diff flux-expr cons-expr)) cons-expr))])
       (not (is-non-negative deriv parameters))) #f]
    
    ;; Otherwise, return true.
    [else #t]))

  (untrace is-real)
  (untrace symbolic-simp)
  (untrace symbolic-simp-rule)
  (untrace symbolic-diff)
  (untrace is-non-negative)

  out)
(trace prove-lax-friedrichs-scalar-1d-local-lipschitz)

;; -------------------------------------------------------------------------
;; Prove hyperbolicity of the Roe (Finite-Volume) Solver for a 1D Scalar PDE
;; -------------------------------------------------------------------------
(define (prove-roe-scalar-1d-hyperbolicity pde
                                           #:nx [nx 200]
                                           #:x0 [x0 0.0]
                                           #:x1 [x1 2.0]
                                           #:t-final [t-final 1.0]
                                           #:cfl [cfl 0.95]
                                           #:init-func [init-func `(cond
                                                                     [(< x 1.0) 1.0]
                                                                     [else 0.0])])
   "Prove that the Roe finite-volume method preserves hyperbolicity for the 1D scalar PDE specified by `pde`. 
  - `nx` : Number of spatial cells.
  - `x0`, `x1` : Domain boundaries.
  - `t-final`: Final time.
  - `cfl`: CFL coefficient.
  - `init-func`: Racket expression for the initial condition, e.g. piecewise constant."

  (define cons-expr (hash-ref pde 'cons-expr))
  (define flux-expr (hash-ref pde 'flux-expr))
  (define parameters (hash-ref pde 'parameters))

  (trace is-real)
  (trace symbolic-simp)
  (trace symbolic-simp-rule)
  (trace symbolic-diff)
  (trace symbolic-roe-function)
  (trace flux-deriv-replace)

  (define flux-deriv (symbolic-simp (symbolic-diff flux-expr cons-expr)))
  
  (define out (cond
    ;; Check whether the CFL coefficient is greater than 0 and less than or equal to 1 (otherwise, return false).
    [(or (<= cfl 0) (> cfl 1)) #f]
    
    ;; Check whether the number of spatial cells is at least 1 and the right domain boundary is set to the right of the left boundary (otherwise, return false)
    [(or (< nx 1) (>= x0 x1)) #f]
    
    ;; Check whether the final simulation time is non-negative (otherwise, return false).
    [(< t-final 0) #f]

    ;; Check whether the simulation parameter(s) correspond to real numbers (otherwise, return false).
    [(not (or (empty? parameters) (andmap (lambda (parameter)
                                            (is-real (list-ref parameter 2) (list cons-expr) parameters)) parameters))) #f]

    ;; Check whether the initial condition(s) correspond to real numbers (otherwise, return false).
    [(not (is-real init-func (list cons-expr) parameters)) #f]
    
    ;; Check whether the Roe function is real (otherwise, return false).
    [(not (is-real (symbolic-roe-function flux-deriv cons-expr) (list
                                                        (string->symbol (string-append (symbol->string cons-expr) "L"))
                                                        (string->symbol (string-append (symbol->string cons-expr) "R"))) parameters)) #f]

    ;; Otherwise, return true.
    [else #t]))

  (untrace is-real)
  (untrace symbolic-simp)
  (untrace symbolic-simp-rule)
  (untrace symbolic-diff)
  (untrace symbolic-roe-function)
  (untrace flux-deriv-replace)
  
  out)
(trace prove-roe-scalar-1d-hyperbolicity)

;; -----------------------------------------------------------------------------------------------
;; Prove flux conservation (jump continuity) of the Roe (Finite-Volume) Solver for a 1D Scalar PDE
;; -----------------------------------------------------------------------------------------------
(define (prove-roe-scalar-1d-flux-conservation pde
                                               #:nx [nx 200]
                                               #:x0 [x0 0.0]
                                               #:x1 [x1 2.0]
                                               #:t-final [t-final 1.0]
                                               #:cfl [cfl 0.95]
                                               #:init-func [init-func `(cond
                                                                         [(< x 1.0) 1.0]
                                                                         [else 0.0])])
   "Prove that the Roe finite-volume method preserves flux conservation (jump continuity) for the 1D scalar PDE specified by `pde`. 
  - `nx` : Number of spatial cells.
  - `x0`, `x1` : Domain boundaries.
  - `t-final`: Final time.
  - `cfl`: CFL coefficient.
  - `init-func`: Racket expression for the initial condition, e.g. piecewise constant."

  (define cons-expr (hash-ref pde 'cons-expr))
  (define flux-expr (hash-ref pde 'flux-expr))
  (define parameters (hash-ref pde 'parameters))

  (trace is-real)
  (trace symbolic-simp)
  (trace symbolic-simp-rule)
  (trace symbolic-diff)
  (trace symbolic-roe-function)
  (trace flux-deriv-replace)

  (define flux-deriv (symbolic-simp (symbolic-diff flux-expr cons-expr)))

  (define roe-jump (symbolic-simp `(* ,(symbolic-roe-function flux-deriv cons-expr) (- ,(string->symbol (string-append (symbol->string cons-expr) "L"))
                                                                                       ,(string->symbol (string-append (symbol->string cons-expr) "R"))))))
  (define flux-jump (symbolic-simp `(- ,(flux-deriv-replace flux-expr cons-expr (string->symbol (string-append (symbol->string cons-expr) "L")))
                                       ,(flux-deriv-replace flux-expr cons-expr (string->symbol (string-append (symbol->string cons-expr) "R"))))))
  
  (define out (cond
    ;; Check whether the CFL coefficient is greater than 0 and less than or equal to 1 (otherwise, return false).
    [(or (<= cfl 0) (> cfl 1)) #f]
    
    ;; Check whether the number of spatial cells is at least 1 and the right domain boundary is set to the right of the left boundary (otherwise, return false)
    [(or (< nx 1) (>= x0 x1)) #f]
    
    ;; Check whether the final simulation time is non-negative (otherwise, return false).
    [(< t-final 0) #f]

    ;; Check whether the simulation parameter(s) correspond to real numbers (otherwise, return false).
    [(not (or (empty? parameters) (andmap (lambda (parameter)
                                            (is-real (list-ref parameter 2) (list cons-expr) parameters)) parameters))) #f]

    ;; Check whether the initial condition(s) correspond to real numbers (otherwise, return false).
    [(not (is-real init-func (list cons-expr) parameters)) #f]
    
    ;; Check whether the jump in the flux function is equal to the product of the Roe function and the jump in the conserved variable (otherwise, return false).
    [(not (equal? roe-jump flux-jump)) #f]

    ;; Otherwise, return true.
    [else #t]))

  (untrace is-real)
  (untrace symbolic-simp)
  (untrace symbolic-simp-rule)
  (untrace symbolic-diff)
  (untrace symbolic-roe-function)
  (untrace flux-deriv-replace)
  
  out)
(trace prove-roe-scalar-1d-flux-conservation)

;; -------------------------------------------------
;; Prove symmetry for a High-Resolution Flux Limiter
;; -------------------------------------------------
(define (prove-flux-limiter-symmetry limiter)
   "Prove that the high-resolution flux limiter specified by `limiter-code` acts symmetrically on forward and backward gradients."

  (define limiter-expr (hash-ref limiter 'limiter-expr))
  (define limiter-ratio (hash-ref limiter 'limiter-ratio))

  (trace variable-transform)
  (trace symbolic-simp)
  (trace symbolic-simp-rule)
  (trace symbolic-simp-positive)
  (trace symbolic-simp-positive-rule)
  
  (define out (cond
    ;; Check whether the symmetry property phi(r) / r = phi(1 / r) holds (otherwise, return false).
    [(not (equal? (symbolic-simp
                   (symbolic-simp-positive (symbolic-simp (symbolic-simp-positive `(/ ,limiter-expr ,limiter-ratio) limiter-ratio)) limiter-ratio))
                  (symbolic-simp
                   (symbolic-simp-positive (symbolic-simp (symbolic-simp-positive (variable-transform limiter-expr limiter-ratio `(/ 1.0 ,limiter-ratio))
                                                                                  limiter-ratio)) limiter-ratio)))) #f]

    ;; Otherwise, return true.
    [else #t]))

  (untrace variable-transform)
  (untrace symbolic-simp)
  (untrace symbolic-simp-rule)
  (untrace symbolic-simp-positive)
  (untrace symbolic-simp-positive-rule)
  
  out)
(trace prove-flux-limiter-symmetry)

;; ---------------------------------------------------------------------------------------
;; Prove second-order TVD (total variation diminishing) for a High-Resolution Flux Limiter
;; ---------------------------------------------------------------------------------------
(define (prove-flux-limiter-tvd limiter)
   "Prove that the high-resolution flux limiter specified by `limiter-code` is second-order TVD (total variation diminishing)."

  (define limiter-expr (hash-ref limiter 'limiter-expr))
  (define limiter-ratio (hash-ref limiter 'limiter-ratio))

  (trace variable-transform)
  (trace symbolic-simp)
  (trace symbolic-simp-rule)
  (trace symbolic-simp-positive)
  (trace symbolic-simp-positive-rule)
  (trace evaluate-limit)
  (trace evaluate-limit-rule)

  (define limiter-convexity (symbolic-simp (symbolic-diff (symbolic-simp-positive (symbolic-simp (symbolic-diff
                                                                                                  (symbolic-simp-positive
                                                                                                   (symbolic-simp limiter-expr) limiter-ratio) limiter-ratio))
                                                                                  limiter-ratio) limiter-ratio)))
  (define limiter-mid (evaluate-limit limiter-expr limiter-ratio 1.0))
  (define limiter-boundary-left (evaluate-limit limiter-expr limiter-ratio 0.0))
  (define limiter-boundary-right (evaluate-limit limiter-expr limiter-ratio 2.0))
  (define limiter-infinity (evaluate-limit limiter-expr limiter-ratio +inf.0))
  
  (define out (cond
    ;; Check whether the limiter function is concave, i.e. that the second derivative of the limiter function is negative (otherwise, return false).
    [(or (not (number? limiter-convexity)) (> limiter-convexity 0.0)) #f]

    ;; Check whether the limiter function limits to 1.0 at the midpoint r = 1.0 (otherwise, return false).
    [(or (not (number? limiter-mid)) (not (equal? limiter-mid 1.0))) #f]

    ;; Check whether the limiter function limits to between 0.0 and 1.0 inclusive at the left (r = 0.0) boundary (otherwise, return false).
    [(or (not (number? limiter-boundary-left)) (> limiter-boundary-left 1.0) (< limiter-boundary-left 0.0)) #f]

    ;; Check whether the limiter function limits to between 1.0 and 2.0 inclusive at the right (r = 2.0) boundary (otherwise, return false).
    [(or (not (number? limiter-boundary-right)) (> limiter-boundary-right 2.0) (< limiter-boundary-right 1.0)) #f]

    ;; Check whether the limiter function limits to less than 2.0 inclusive as r approaches +infinity (otherwise, return false).
    [(or (not (number? limiter-infinity)) (> limiter-infinity 2.0)) #f]

    ;; Otherwise, return true.
    [else #t]))

  (untrace variable-transform)
  (untrace symbolic-simp)
  (untrace symbolic-simp-rule)
  (untrace symbolic-simp-positive)
  (untrace symbolic-simp-positive-rule)
  (untrace evaluate-limit)
  (untrace evaluate-limit-rule)
  
  out)
(trace prove-flux-limiter-tvd)