#lang racket

(require "../prover_core.rkt")
(require "../prover_vector.rkt")

 (prove-lax-friedrichs-vector2-1d-local-lipschitz '#hash((cons-exprs . (Ey Bz)) (flux-exprs . ((* (* c c) Bz) Ey)) (max-speed-exprs . ((abs c) (abs c))) (name . "maxwell_EyBz") (parameters . ((define c 1.0) (define e_fact 1.0) (define b_fact 1.0)))) #:cfl 0.95 #:init-funcs '(0.0 (cond ((< x 0.0) 0.5) (else -0.5))) #:nx 200 #:t-final 1.0 #:x0 -1.5 #:x1 1.5)
  (symbolic-hessian '(* (* c c) Bz) '(Ey Bz))
   (symbolic-gradient '(* (* c c) Bz) '(Ey Bz))
    (symbolic-diff '(* (* c c) Bz) 'Ey)
     (symbolic-diff '(* c c) 'Ey)
      (symbolic-diff 'c 'Ey)
      0.0
      (symbolic-diff 'c 'Ey)
      0.0
     '(+ (* 0.0 c) (* c 0.0))
     (symbolic-diff 'Bz 'Ey)
     0.0
    '(+ (* (+ (* 0.0 c) (* c 0.0)) Bz) (* (* c c) 0.0))
    (symbolic-simp '(+ (* (+ (* 0.0 c) (* c 0.0)) Bz) (* (* c c) 0.0)))
     (symbolic-simp-rule '(+ (* (+ (* 0.0 c) (* c 0.0)) Bz) (* (* c c) 0.0)))
      (symbolic-simp-rule '(* (+ (* 0.0 c) (* c 0.0)) Bz))
       (symbolic-simp-rule '(+ (* 0.0 c) (* c 0.0)))
        (symbolic-simp-rule '(* 0.0 c))
        0.0
        (symbolic-simp-rule '(* c 0.0))
        '(* 0.0 c)
       '(+ 0.0 (* 0.0 c))
       (symbolic-simp-rule 'Bz)
       'Bz
      '(* (+ 0.0 (* 0.0 c)) Bz)
      (symbolic-simp-rule '(* (* c c) 0.0))
      '(* c (* c 0.0))
     '(+ (* (+ 0.0 (* 0.0 c)) Bz) (* c (* c 0.0)))
    (symbolic-simp '(+ (* (+ 0.0 (* 0.0 c)) Bz) (* c (* c 0.0))))
     (symbolic-simp-rule '(+ (* (+ 0.0 (* 0.0 c)) Bz) (* c (* c 0.0))))
      (symbolic-simp-rule '(* (+ 0.0 (* 0.0 c)) Bz))
       (symbolic-simp-rule '(+ 0.0 (* 0.0 c)))
       '(* 0.0 c)
       (symbolic-simp-rule 'Bz)
       'Bz
      '(* (* 0.0 c) Bz)
      (symbolic-simp-rule '(* c (* c 0.0)))
       (symbolic-simp-rule 'c)
       'c
       (symbolic-simp-rule '(* c 0.0))
       '(* 0.0 c)
      '(* c (* 0.0 c))
     '(+ (* (* 0.0 c) Bz) (* c (* 0.0 c)))
    (symbolic-simp '(+ (* (* 0.0 c) Bz) (* c (* 0.0 c))))
     (symbolic-simp-rule '(+ (* (* 0.0 c) Bz) (* c (* 0.0 c))))
      (symbolic-simp-rule '(* (* 0.0 c) Bz))
      '(* 0.0 (* c Bz))
      (symbolic-simp-rule '(* c (* 0.0 c)))
      '(* 0.0 (* c c))
     '(+ (* 0.0 (* c Bz)) (* 0.0 (* c c)))
    (symbolic-simp '(+ (* 0.0 (* c Bz)) (* 0.0 (* c c))))
     (symbolic-simp-rule '(+ (* 0.0 (* c Bz)) (* 0.0 (* c c))))
      (symbolic-simp-rule '(* 0.0 (* c Bz)))
      0.0
      (symbolic-simp-rule '(* 0.0 (* c c)))
      0.0
     '(+ 0.0 0.0)
    (symbolic-simp '(+ 0.0 0.0))
     (symbolic-simp-rule '(+ 0.0 0.0))
     0.0
    (symbolic-simp 0.0)
     (symbolic-simp-rule 0.0)
     0.0
    0.0
    (symbolic-diff '(* (* c c) Bz) 'Bz)
     (symbolic-diff '(* c c) 'Bz)
      (symbolic-diff 'c 'Bz)
      0.0
      (symbolic-diff 'c 'Bz)
      0.0
     '(+ (* 0.0 c) (* c 0.0))
     (symbolic-diff 'Bz 'Bz)
     1.0
    '(+ (* (+ (* 0.0 c) (* c 0.0)) Bz) (* (* c c) 1.0))
    (symbolic-simp '(+ (* (+ (* 0.0 c) (* c 0.0)) Bz) (* (* c c) 1.0)))
     (symbolic-simp-rule '(+ (* (+ (* 0.0 c) (* c 0.0)) Bz) (* (* c c) 1.0)))
      (symbolic-simp-rule '(* (+ (* 0.0 c) (* c 0.0)) Bz))
       (symbolic-simp-rule '(+ (* 0.0 c) (* c 0.0)))
        (symbolic-simp-rule '(* 0.0 c))
        0.0
        (symbolic-simp-rule '(* c 0.0))
        '(* 0.0 c)
       '(+ 0.0 (* 0.0 c))
       (symbolic-simp-rule 'Bz)
       'Bz
      '(* (+ 0.0 (* 0.0 c)) Bz)
      (symbolic-simp-rule '(* (* c c) 1.0))
      '(* c (* c 1.0))
     '(+ (* (+ 0.0 (* 0.0 c)) Bz) (* c (* c 1.0)))
    (symbolic-simp '(+ (* (+ 0.0 (* 0.0 c)) Bz) (* c (* c 1.0))))
     (symbolic-simp-rule '(+ (* (+ 0.0 (* 0.0 c)) Bz) (* c (* c 1.0))))
      (symbolic-simp-rule '(* (+ 0.0 (* 0.0 c)) Bz))
       (symbolic-simp-rule '(+ 0.0 (* 0.0 c)))
       '(* 0.0 c)
       (symbolic-simp-rule 'Bz)
       'Bz
      '(* (* 0.0 c) Bz)
      (symbolic-simp-rule '(* c (* c 1.0)))
       (symbolic-simp-rule 'c)
       'c
       (symbolic-simp-rule '(* c 1.0))
       '(* 1.0 c)
      '(* c (* 1.0 c))
     '(+ (* (* 0.0 c) Bz) (* c (* 1.0 c)))
    (symbolic-simp '(+ (* (* 0.0 c) Bz) (* c (* 1.0 c))))
     (symbolic-simp-rule '(+ (* (* 0.0 c) Bz) (* c (* 1.0 c))))
      (symbolic-simp-rule '(* (* 0.0 c) Bz))
      '(* 0.0 (* c Bz))
      (symbolic-simp-rule '(* c (* 1.0 c)))
      '(* 1.0 (* c c))
     '(+ (* 0.0 (* c Bz)) (* 1.0 (* c c)))
    (symbolic-simp '(+ (* 0.0 (* c Bz)) (* 1.0 (* c c))))
     (symbolic-simp-rule '(+ (* 0.0 (* c Bz)) (* 1.0 (* c c))))
      (symbolic-simp-rule '(* 0.0 (* c Bz)))
      0.0
      (symbolic-simp-rule '(* 1.0 (* c c)))
      '(* c c)
     '(+ 0.0 (* c c))
    (symbolic-simp '(+ 0.0 (* c c)))
     (symbolic-simp-rule '(+ 0.0 (* c c)))
     '(* c c)
    (symbolic-simp '(* c c))
     (symbolic-simp-rule '(* c c))
      (symbolic-simp-rule 'c)
      'c
      (symbolic-simp-rule 'c)
      'c
     '(* c c)
    '(* c c)
   '(0.0 (* c c))
  (symbolic-jacobian '(0.0 (* c c)) '(Ey Bz))
   (symbolic-diff 0.0 'Ey)
   0.0
   (symbolic-simp 0.0)
    (symbolic-simp-rule 0.0)
    0.0
   0.0
   (symbolic-diff 0.0 'Bz)
   0.0
   (symbolic-simp 0.0)
    (symbolic-simp-rule 0.0)
    0.0
   0.0
   (symbolic-diff '(* c c) 'Ey)
    (symbolic-diff 'c 'Ey)
    0.0
    (symbolic-diff 'c 'Ey)
    0.0
   '(+ (* 0.0 c) (* c 0.0))
   (symbolic-simp '(+ (* 0.0 c) (* c 0.0)))
    (symbolic-simp-rule '(+ (* 0.0 c) (* c 0.0)))
     (symbolic-simp-rule '(* 0.0 c))
     0.0
     (symbolic-simp-rule '(* c 0.0))
     '(* 0.0 c)
    '(+ 0.0 (* 0.0 c))
   (symbolic-simp '(+ 0.0 (* 0.0 c)))
    (symbolic-simp-rule '(+ 0.0 (* 0.0 c)))
    '(* 0.0 c)
   (symbolic-simp '(* 0.0 c))
    (symbolic-simp-rule '(* 0.0 c))
    0.0
   (symbolic-simp 0.0)
    (symbolic-simp-rule 0.0)
    0.0
   0.0
   (symbolic-diff '(* c c) 'Bz)
    (symbolic-diff 'c 'Bz)
    0.0
    (symbolic-diff 'c 'Bz)
    0.0
   '(+ (* 0.0 c) (* c 0.0))
   (symbolic-simp '(+ (* 0.0 c) (* c 0.0)))
    (symbolic-simp-rule '(+ (* 0.0 c) (* c 0.0)))
     (symbolic-simp-rule '(* 0.0 c))
     0.0
     (symbolic-simp-rule '(* c 0.0))
     '(* 0.0 c)
    '(+ 0.0 (* 0.0 c))
   (symbolic-simp '(+ 0.0 (* 0.0 c)))
    (symbolic-simp-rule '(+ 0.0 (* 0.0 c)))
    '(* 0.0 c)
   (symbolic-simp '(* 0.0 c))
    (symbolic-simp-rule '(* 0.0 c))
    0.0
   (symbolic-simp 0.0)
    (symbolic-simp-rule 0.0)
    0.0
   0.0
  '((0.0 0.0) (0.0 0.0))
  (symbolic-hessian 'Ey '(Ey Bz))
   (symbolic-gradient 'Ey '(Ey Bz))
    (symbolic-diff 'Ey 'Ey)
    1.0
    (symbolic-simp 1.0)
     (symbolic-simp-rule 1.0)
     1.0
    1.0
    (symbolic-diff 'Ey 'Bz)
    0.0
    (symbolic-simp 0.0)
     (symbolic-simp-rule 0.0)
     0.0
    0.0
   '(1.0 0.0)
  (symbolic-jacobian '(1.0 0.0) '(Ey Bz))
   (symbolic-diff 1.0 'Ey)
   0.0
   (symbolic-simp 0.0)
    (symbolic-simp-rule 0.0)
    0.0
   0.0
   (symbolic-diff 1.0 'Bz)
   0.0
   (symbolic-simp 0.0)
    (symbolic-simp-rule 0.0)
    0.0
   0.0
   (symbolic-diff 0.0 'Ey)
   0.0
   (symbolic-simp 0.0)
    (symbolic-simp-rule 0.0)
    0.0
   0.0
   (symbolic-diff 0.0 'Bz)
   0.0
   (symbolic-simp 0.0)
    (symbolic-simp-rule 0.0)
    0.0
   0.0
  '((0.0 0.0) (0.0 0.0))
  (symbolic-eigvals2 '((0.0 0.0) (0.0 0.0)))
  '(0.0 0.0)
  (symbolic-eigvals2 '((0.0 0.0) (0.0 0.0)))
  '(0.0 0.0)
  (symbolic-simp 0.0)
   (symbolic-simp-rule 0.0)
   0.0
  0.0
  (symbolic-simp 0.0)
   (symbolic-simp-rule 0.0)
   0.0
  0.0
  (symbolic-simp 0.0)
   (symbolic-simp-rule 0.0)
   0.0
  0.0
  (symbolic-simp 0.0)
   (symbolic-simp-rule 0.0)
   0.0
  0.0
  (is-real 1.0 '(Ey Bz) '((define c 1.0) (define e_fact 1.0) (define b_fact 1.0)))
  #t
  (is-real 1.0 '(Ey Bz) '((define c 1.0) (define e_fact 1.0) (define b_fact 1.0)))
  #t
  (is-real 1.0 '(Ey Bz) '((define c 1.0) (define e_fact 1.0) (define b_fact 1.0)))
  #t
  (is-real 0.0 '(Ey Bz) '((define c 1.0) (define e_fact 1.0) (define b_fact 1.0)))
  #t
  (is-real '(cond ((< x 0.0) 0.5) (else -0.5)) '(Ey Bz) '((define c 1.0) (define e_fact 1.0) (define b_fact 1.0)))
   (is-real 0.5 '(Ey Bz) '((define c 1.0) (define e_fact 1.0) (define b_fact 1.0)))
   #t
  (is-real -0.5 '(Ey Bz) '((define c 1.0) (define e_fact 1.0) (define b_fact 1.0)))
  #t
  (is-non-negative 0.0 '((define c 1.0) (define e_fact 1.0) (define b_fact 1.0)))
  #t
  (is-non-negative 0.0 '((define c 1.0) (define e_fact 1.0) (define b_fact 1.0)))
  #t
  (is-non-negative 0.0 '((define c 1.0) (define e_fact 1.0) (define b_fact 1.0)))
  #t
  (is-non-negative 0.0 '((define c 1.0) (define e_fact 1.0) (define b_fact 1.0)))
  #t
 #t
