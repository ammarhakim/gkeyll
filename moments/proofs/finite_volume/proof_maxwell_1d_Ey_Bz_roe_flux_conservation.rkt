#lang racket

(require "../prover_core.rkt")
(require "../prover_vector.rkt")

 (prove-roe-vector2-1d-flux-conservation '#hash((cons-exprs . (Ey Bz)) (flux-exprs . ((* (* c c) Bz) Ey)) (max-speed-exprs . ((abs c) (abs c))) (name . "maxwell_EyBz") (parameters . ((define c 1.0) (define e_fact 1.0) (define b_fact 1.0)))) #:cfl 0.95 #:init-funcs '(0.0 (cond ((< x 0.0) 0.5) (else -0.5))) #:nx 200 #:t-final 1.0 #:x0 -1.5 #:x1 1.5)
  (symbolic-jacobian '((* (* c c) Bz) Ey) '(Ey Bz))
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
  '((0.0 (* c c)) (1.0 0.0))
  (symbolic-roe-matrix '((0.0 (* c c)) (1.0 0.0)) '(Ey Bz))
   (flux-deriv-replace 0.0 'Ey 'EyL)
   0.0
   (flux-deriv-replace 0.0 'Bz 'BzL)
   0.0
   (flux-deriv-replace 0.0 'Ey 'EyR)
   0.0
   (flux-deriv-replace 0.0 'Bz 'BzR)
   0.0
   (symbolic-simp '(+ (* 0.5 0.0) (* 0.5 0.0)))
    (symbolic-simp-rule '(+ (* 0.5 0.0) (* 0.5 0.0)))
    '(* (+ 0.5 0.5) 0.0)
   (symbolic-simp '(* (+ 0.5 0.5) 0.0))
    (symbolic-simp-rule '(* (+ 0.5 0.5) 0.0))
    '(* 0.0 (+ 0.5 0.5))
   (symbolic-simp '(* 0.0 (+ 0.5 0.5)))
    (symbolic-simp-rule '(* 0.0 (+ 0.5 0.5)))
    0.0
   (symbolic-simp 0.0)
    (symbolic-simp-rule 0.0)
    0.0
   0.0
   (flux-deriv-replace '(* c c) 'Ey 'EyL)
    (flux-deriv-replace 'c 'Ey 'EyL)
    'c
    (flux-deriv-replace 'c 'Ey 'EyL)
    'c
   '(* c c)
   (flux-deriv-replace '(* c c) 'Bz 'BzL)
    (flux-deriv-replace 'c 'Bz 'BzL)
    'c
    (flux-deriv-replace 'c 'Bz 'BzL)
    'c
   '(* c c)
   (flux-deriv-replace '(* c c) 'Ey 'EyR)
    (flux-deriv-replace 'c 'Ey 'EyR)
    'c
    (flux-deriv-replace 'c 'Ey 'EyR)
    'c
   '(* c c)
   (flux-deriv-replace '(* c c) 'Bz 'BzR)
    (flux-deriv-replace 'c 'Bz 'BzR)
    'c
    (flux-deriv-replace 'c 'Bz 'BzR)
    'c
   '(* c c)
   (symbolic-simp '(+ (* 0.5 (* c c)) (* 0.5 (* c c))))
    (symbolic-simp-rule '(+ (* 0.5 (* c c)) (* 0.5 (* c c))))
    '(* (+ 0.5 0.5) (* c c))
   (symbolic-simp '(* (+ 0.5 0.5) (* c c)))
    (symbolic-simp-rule '(* (+ 0.5 0.5) (* c c)))
     (symbolic-simp-rule '(+ 0.5 0.5))
     1.0
     (symbolic-simp-rule '(* c c))
      (symbolic-simp-rule 'c)
      'c
      (symbolic-simp-rule 'c)
      'c
     '(* c c)
    '(* 1.0 (* c c))
   (symbolic-simp '(* 1.0 (* c c)))
    (symbolic-simp-rule '(* 1.0 (* c c)))
    '(* c c)
   (symbolic-simp '(* c c))
    (symbolic-simp-rule '(* c c))
     (symbolic-simp-rule 'c)
     'c
     (symbolic-simp-rule 'c)
     'c
    '(* c c)
   '(* c c)
   (flux-deriv-replace 1.0 'Ey 'EyL)
   1.0
   (flux-deriv-replace 1.0 'Bz 'BzL)
   1.0
   (flux-deriv-replace 1.0 'Ey 'EyR)
   1.0
   (flux-deriv-replace 1.0 'Bz 'BzR)
   1.0
   (symbolic-simp '(+ (* 0.5 1.0) (* 0.5 1.0)))
    (symbolic-simp-rule '(+ (* 0.5 1.0) (* 0.5 1.0)))
    '(* (+ 0.5 0.5) 1.0)
   (symbolic-simp '(* (+ 0.5 0.5) 1.0))
    (symbolic-simp-rule '(* (+ 0.5 0.5) 1.0))
    '(* 1.0 (+ 0.5 0.5))
   (symbolic-simp '(* 1.0 (+ 0.5 0.5)))
    (symbolic-simp-rule '(* 1.0 (+ 0.5 0.5)))
    '(+ 0.5 0.5)
   (symbolic-simp '(+ 0.5 0.5))
    (symbolic-simp-rule '(+ 0.5 0.5))
    1.0
   (symbolic-simp 1.0)
    (symbolic-simp-rule 1.0)
    1.0
   1.0
   (flux-deriv-replace 0.0 'Ey 'EyL)
   0.0
   (flux-deriv-replace 0.0 'Bz 'BzL)
   0.0
   (flux-deriv-replace 0.0 'Ey 'EyR)
   0.0
   (flux-deriv-replace 0.0 'Bz 'BzR)
   0.0
   (symbolic-simp '(+ (* 0.5 0.0) (* 0.5 0.0)))
    (symbolic-simp-rule '(+ (* 0.5 0.0) (* 0.5 0.0)))
    '(* (+ 0.5 0.5) 0.0)
   (symbolic-simp '(* (+ 0.5 0.5) 0.0))
    (symbolic-simp-rule '(* (+ 0.5 0.5) 0.0))
    '(* 0.0 (+ 0.5 0.5))
   (symbolic-simp '(* 0.0 (+ 0.5 0.5)))
    (symbolic-simp-rule '(* 0.0 (+ 0.5 0.5)))
    0.0
   (symbolic-simp 0.0)
    (symbolic-simp-rule 0.0)
    0.0
   0.0
  '((0.0 (* c c)) (1.0 0.0))
  (symbolic-simp '(- EyL EyR))
   (symbolic-simp-rule '(- EyL EyR))
    (symbolic-simp-rule 'EyL)
    'EyL
    (symbolic-simp-rule 'EyR)
    'EyR
   '(- EyL EyR)
  '(- EyL EyR)
  (symbolic-simp '(- BzL BzR))
   (symbolic-simp-rule '(- BzL BzR))
    (symbolic-simp-rule 'BzL)
    'BzL
    (symbolic-simp-rule 'BzR)
    'BzR
   '(- BzL BzR)
  '(- BzL BzR)
  (symbolic-simp '(+ (* 0.0 (- EyL EyR)) (* (* c c) (- BzL BzR))))
   (symbolic-simp-rule '(+ (* 0.0 (- EyL EyR)) (* (* c c) (- BzL BzR))))
    (symbolic-simp-rule '(* 0.0 (- EyL EyR)))
    0.0
    (symbolic-simp-rule '(* (* c c) (- BzL BzR)))
    '(* c (* c (- BzL BzR)))
   '(+ 0.0 (* c (* c (- BzL BzR))))
  (symbolic-simp '(+ 0.0 (* c (* c (- BzL BzR)))))
   (symbolic-simp-rule '(+ 0.0 (* c (* c (- BzL BzR)))))
   '(* c (* c (- BzL BzR)))
  (symbolic-simp '(* c (* c (- BzL BzR))))
   (symbolic-simp-rule '(* c (* c (- BzL BzR))))
    (symbolic-simp-rule 'c)
    'c
    (symbolic-simp-rule '(* c (- BzL BzR)))
     (symbolic-simp-rule 'c)
     'c
     (symbolic-simp-rule '(- BzL BzR))
      (symbolic-simp-rule 'BzL)
      'BzL
      (symbolic-simp-rule 'BzR)
      'BzR
     '(- BzL BzR)
    '(* c (- BzL BzR))
   '(* c (* c (- BzL BzR)))
  '(* c (* c (- BzL BzR)))
  (symbolic-simp '(+ (* 1.0 (- EyL EyR)) (* 0.0 (- BzL BzR))))
   (symbolic-simp-rule '(+ (* 1.0 (- EyL EyR)) (* 0.0 (- BzL BzR))))
    (symbolic-simp-rule '(* 1.0 (- EyL EyR)))
    '(- EyL EyR)
    (symbolic-simp-rule '(* 0.0 (- BzL BzR)))
    0.0
   '(+ (- EyL EyR) 0.0)
  (symbolic-simp '(+ (- EyL EyR) 0.0))
   (symbolic-simp-rule '(+ (- EyL EyR) 0.0))
   '(+ 0.0 (- EyL EyR))
  (symbolic-simp '(+ 0.0 (- EyL EyR)))
   (symbolic-simp-rule '(+ 0.0 (- EyL EyR)))
   '(- EyL EyR)
  (symbolic-simp '(- EyL EyR))
   (symbolic-simp-rule '(- EyL EyR))
    (symbolic-simp-rule 'EyL)
    'EyL
    (symbolic-simp-rule 'EyR)
    'EyR
   '(- EyL EyR)
  '(- EyL EyR)
  (flux-deriv-replace '(* (* c c) Bz) 'Ey 'EyL)
   (flux-deriv-replace '(* c c) 'Ey 'EyL)
    (flux-deriv-replace 'c 'Ey 'EyL)
    'c
    (flux-deriv-replace 'c 'Ey 'EyL)
    'c
   '(* c c)
   (flux-deriv-replace 'Bz 'Ey 'EyL)
   'Bz
  '(* (* c c) Bz)
  (flux-deriv-replace '(* (* c c) Bz) 'Bz 'BzL)
   (flux-deriv-replace '(* c c) 'Bz 'BzL)
    (flux-deriv-replace 'c 'Bz 'BzL)
    'c
    (flux-deriv-replace 'c 'Bz 'BzL)
    'c
   '(* c c)
   (flux-deriv-replace 'Bz 'Bz 'BzL)
   'BzL
  '(* (* c c) BzL)
  (flux-deriv-replace '(* (* c c) Bz) 'Ey 'EyR)
   (flux-deriv-replace '(* c c) 'Ey 'EyR)
    (flux-deriv-replace 'c 'Ey 'EyR)
    'c
    (flux-deriv-replace 'c 'Ey 'EyR)
    'c
   '(* c c)
   (flux-deriv-replace 'Bz 'Ey 'EyR)
   'Bz
  '(* (* c c) Bz)
  (flux-deriv-replace '(* (* c c) Bz) 'Bz 'BzR)
   (flux-deriv-replace '(* c c) 'Bz 'BzR)
    (flux-deriv-replace 'c 'Bz 'BzR)
    'c
    (flux-deriv-replace 'c 'Bz 'BzR)
    'c
   '(* c c)
   (flux-deriv-replace 'Bz 'Bz 'BzR)
   'BzR
  '(* (* c c) BzR)
  (symbolic-simp '(- (* (* c c) BzL) (* (* c c) BzR)))
   (symbolic-simp-rule '(- (* (* c c) BzL) (* (* c c) BzR)))
   '(* (* c c) (- BzL BzR))
  (symbolic-simp '(* (* c c) (- BzL BzR)))
   (symbolic-simp-rule '(* (* c c) (- BzL BzR)))
   '(* c (* c (- BzL BzR)))
  (symbolic-simp '(* c (* c (- BzL BzR))))
   (symbolic-simp-rule '(* c (* c (- BzL BzR))))
    (symbolic-simp-rule 'c)
    'c
    (symbolic-simp-rule '(* c (- BzL BzR)))
     (symbolic-simp-rule 'c)
     'c
     (symbolic-simp-rule '(- BzL BzR))
      (symbolic-simp-rule 'BzL)
      'BzL
      (symbolic-simp-rule 'BzR)
      'BzR
     '(- BzL BzR)
    '(* c (- BzL BzR))
   '(* c (* c (- BzL BzR)))
  '(* c (* c (- BzL BzR)))
  (flux-deriv-replace 'Ey 'Ey 'EyL)
  'EyL
  (flux-deriv-replace 'EyL 'Bz 'BzL)
  'EyL
  (flux-deriv-replace 'Ey 'Ey 'EyR)
  'EyR
  (flux-deriv-replace 'EyR 'Bz 'BzR)
  'EyR
  (symbolic-simp '(- EyL EyR))
   (symbolic-simp-rule '(- EyL EyR))
    (symbolic-simp-rule 'EyL)
    'EyL
    (symbolic-simp-rule 'EyR)
    'EyR
   '(- EyL EyR)
  '(- EyL EyR)
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
 #t
