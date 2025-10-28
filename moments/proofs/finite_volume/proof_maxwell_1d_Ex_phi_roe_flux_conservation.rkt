#lang racket

(require "../prover_core.rkt")
(require "../prover_vector.rkt")

 (prove-roe-vector2-1d-flux-conservation '#hash((cons-exprs . (Ex phi)) (flux-exprs . ((* e_fact (* (* c c) phi)) (* e_fact Ex))) (max-speed-exprs . ((abs (* c e_fact)) (abs (* c e_fact)))) (name . "maxwell_Exphi") (parameters . ((define c 1.0) (define e_fact 1.0) (define b_fact 1.0)))) #:cfl 0.95 #:init-funcs '(0.0 (cond ((< x 0.0) 0.5) (else -0.5))) #:nx 200 #:t-final 1.0 #:x0 -1.5 #:x1 1.5)
  (symbolic-jacobian '((* e_fact (* (* c c) phi)) (* e_fact Ex)) '(Ex phi))
   (symbolic-diff '(* e_fact (* (* c c) phi)) 'Ex)
    (symbolic-diff 'e_fact 'Ex)
    0.0
    (symbolic-diff '(* (* c c) phi) 'Ex)
     (symbolic-diff '(* c c) 'Ex)
      (symbolic-diff 'c 'Ex)
      0.0
      (symbolic-diff 'c 'Ex)
      0.0
     '(+ (* 0.0 c) (* c 0.0))
     (symbolic-diff 'phi 'Ex)
     0.0
    '(+ (* (+ (* 0.0 c) (* c 0.0)) phi) (* (* c c) 0.0))
   '(+ (* 0.0 (* (* c c) phi)) (* e_fact (+ (* (+ (* 0.0 c) (* c 0.0)) phi) (* (* c c) 0.0))))
   (symbolic-simp '(+ (* 0.0 (* (* c c) phi)) (* e_fact (+ (* (+ (* 0.0 c) (* c 0.0)) phi) (* (* c c) 0.0)))))
    (symbolic-simp-rule '(+ (* 0.0 (* (* c c) phi)) (* e_fact (+ (* (+ (* 0.0 c) (* c 0.0)) phi) (* (* c c) 0.0)))))
     (symbolic-simp-rule '(* 0.0 (* (* c c) phi)))
     0.0
     (symbolic-simp-rule '(* e_fact (+ (* (+ (* 0.0 c) (* c 0.0)) phi) (* (* c c) 0.0))))
     '(+ (* e_fact (* (+ (* 0.0 c) (* c 0.0)) phi)) (* e_fact (* (* c c) 0.0)))
    '(+ 0.0 (+ (* e_fact (* (+ (* 0.0 c) (* c 0.0)) phi)) (* e_fact (* (* c c) 0.0))))
   (symbolic-simp '(+ 0.0 (+ (* e_fact (* (+ (* 0.0 c) (* c 0.0)) phi)) (* e_fact (* (* c c) 0.0)))))
    (symbolic-simp-rule '(+ 0.0 (+ (* e_fact (* (+ (* 0.0 c) (* c 0.0)) phi)) (* e_fact (* (* c c) 0.0)))))
    '(+ (* e_fact (* (+ (* 0.0 c) (* c 0.0)) phi)) (* e_fact (* (* c c) 0.0)))
   (symbolic-simp '(+ (* e_fact (* (+ (* 0.0 c) (* c 0.0)) phi)) (* e_fact (* (* c c) 0.0))))
    (symbolic-simp-rule '(+ (* e_fact (* (+ (* 0.0 c) (* c 0.0)) phi)) (* e_fact (* (* c c) 0.0))))
     (symbolic-simp-rule '(* e_fact (* (+ (* 0.0 c) (* c 0.0)) phi)))
      (symbolic-simp-rule 'e_fact)
      'e_fact
      (symbolic-simp-rule '(* (+ (* 0.0 c) (* c 0.0)) phi))
       (symbolic-simp-rule '(+ (* 0.0 c) (* c 0.0)))
        (symbolic-simp-rule '(* 0.0 c))
        0.0
        (symbolic-simp-rule '(* c 0.0))
        '(* 0.0 c)
       '(+ 0.0 (* 0.0 c))
       (symbolic-simp-rule 'phi)
       'phi
      '(* (+ 0.0 (* 0.0 c)) phi)
     '(* e_fact (* (+ 0.0 (* 0.0 c)) phi))
     (symbolic-simp-rule '(* e_fact (* (* c c) 0.0)))
      (symbolic-simp-rule 'e_fact)
      'e_fact
      (symbolic-simp-rule '(* (* c c) 0.0))
      '(* c (* c 0.0))
     '(* e_fact (* c (* c 0.0)))
    '(+ (* e_fact (* (+ 0.0 (* 0.0 c)) phi)) (* e_fact (* c (* c 0.0))))
   (symbolic-simp '(+ (* e_fact (* (+ 0.0 (* 0.0 c)) phi)) (* e_fact (* c (* c 0.0)))))
    (symbolic-simp-rule '(+ (* e_fact (* (+ 0.0 (* 0.0 c)) phi)) (* e_fact (* c (* c 0.0)))))
     (symbolic-simp-rule '(* e_fact (* (+ 0.0 (* 0.0 c)) phi)))
      (symbolic-simp-rule 'e_fact)
      'e_fact
      (symbolic-simp-rule '(* (+ 0.0 (* 0.0 c)) phi))
       (symbolic-simp-rule '(+ 0.0 (* 0.0 c)))
       '(* 0.0 c)
       (symbolic-simp-rule 'phi)
       'phi
      '(* (* 0.0 c) phi)
     '(* e_fact (* (* 0.0 c) phi))
     (symbolic-simp-rule '(* e_fact (* c (* c 0.0))))
      (symbolic-simp-rule 'e_fact)
      'e_fact
      (symbolic-simp-rule '(* c (* c 0.0)))
       (symbolic-simp-rule 'c)
       'c
       (symbolic-simp-rule '(* c 0.0))
       '(* 0.0 c)
      '(* c (* 0.0 c))
     '(* e_fact (* c (* 0.0 c)))
    '(+ (* e_fact (* (* 0.0 c) phi)) (* e_fact (* c (* 0.0 c))))
   (symbolic-simp '(+ (* e_fact (* (* 0.0 c) phi)) (* e_fact (* c (* 0.0 c)))))
    (symbolic-simp-rule '(+ (* e_fact (* (* 0.0 c) phi)) (* e_fact (* c (* 0.0 c)))))
     (symbolic-simp-rule '(* e_fact (* (* 0.0 c) phi)))
      (symbolic-simp-rule 'e_fact)
      'e_fact
      (symbolic-simp-rule '(* (* 0.0 c) phi))
      '(* 0.0 (* c phi))
     '(* e_fact (* 0.0 (* c phi)))
     (symbolic-simp-rule '(* e_fact (* c (* 0.0 c))))
      (symbolic-simp-rule 'e_fact)
      'e_fact
      (symbolic-simp-rule '(* c (* 0.0 c)))
      '(* 0.0 (* c c))
     '(* e_fact (* 0.0 (* c c)))
    '(+ (* e_fact (* 0.0 (* c phi))) (* e_fact (* 0.0 (* c c))))
   (symbolic-simp '(+ (* e_fact (* 0.0 (* c phi))) (* e_fact (* 0.0 (* c c)))))
    (symbolic-simp-rule '(+ (* e_fact (* 0.0 (* c phi))) (* e_fact (* 0.0 (* c c)))))
     (symbolic-simp-rule '(* e_fact (* 0.0 (* c phi))))
     '(* 0.0 (* e_fact (* c phi)))
     (symbolic-simp-rule '(* e_fact (* 0.0 (* c c))))
     '(* 0.0 (* e_fact (* c c)))
    '(+ (* 0.0 (* e_fact (* c phi))) (* 0.0 (* e_fact (* c c))))
   (symbolic-simp '(+ (* 0.0 (* e_fact (* c phi))) (* 0.0 (* e_fact (* c c)))))
    (symbolic-simp-rule '(+ (* 0.0 (* e_fact (* c phi))) (* 0.0 (* e_fact (* c c)))))
     (symbolic-simp-rule '(* 0.0 (* e_fact (* c phi))))
     0.0
     (symbolic-simp-rule '(* 0.0 (* e_fact (* c c))))
     0.0
    '(+ 0.0 0.0)
   (symbolic-simp '(+ 0.0 0.0))
    (symbolic-simp-rule '(+ 0.0 0.0))
    0.0
   (symbolic-simp 0.0)
    (symbolic-simp-rule 0.0)
    0.0
   0.0
   (symbolic-diff '(* e_fact (* (* c c) phi)) 'phi)
    (symbolic-diff 'e_fact 'phi)
    0.0
    (symbolic-diff '(* (* c c) phi) 'phi)
     (symbolic-diff '(* c c) 'phi)
      (symbolic-diff 'c 'phi)
      0.0
      (symbolic-diff 'c 'phi)
      0.0
     '(+ (* 0.0 c) (* c 0.0))
     (symbolic-diff 'phi 'phi)
     1.0
    '(+ (* (+ (* 0.0 c) (* c 0.0)) phi) (* (* c c) 1.0))
   '(+ (* 0.0 (* (* c c) phi)) (* e_fact (+ (* (+ (* 0.0 c) (* c 0.0)) phi) (* (* c c) 1.0))))
   (symbolic-simp '(+ (* 0.0 (* (* c c) phi)) (* e_fact (+ (* (+ (* 0.0 c) (* c 0.0)) phi) (* (* c c) 1.0)))))
    (symbolic-simp-rule '(+ (* 0.0 (* (* c c) phi)) (* e_fact (+ (* (+ (* 0.0 c) (* c 0.0)) phi) (* (* c c) 1.0)))))
     (symbolic-simp-rule '(* 0.0 (* (* c c) phi)))
     0.0
     (symbolic-simp-rule '(* e_fact (+ (* (+ (* 0.0 c) (* c 0.0)) phi) (* (* c c) 1.0))))
     '(+ (* e_fact (* (+ (* 0.0 c) (* c 0.0)) phi)) (* e_fact (* (* c c) 1.0)))
    '(+ 0.0 (+ (* e_fact (* (+ (* 0.0 c) (* c 0.0)) phi)) (* e_fact (* (* c c) 1.0))))
   (symbolic-simp '(+ 0.0 (+ (* e_fact (* (+ (* 0.0 c) (* c 0.0)) phi)) (* e_fact (* (* c c) 1.0)))))
    (symbolic-simp-rule '(+ 0.0 (+ (* e_fact (* (+ (* 0.0 c) (* c 0.0)) phi)) (* e_fact (* (* c c) 1.0)))))
    '(+ (* e_fact (* (+ (* 0.0 c) (* c 0.0)) phi)) (* e_fact (* (* c c) 1.0)))
   (symbolic-simp '(+ (* e_fact (* (+ (* 0.0 c) (* c 0.0)) phi)) (* e_fact (* (* c c) 1.0))))
    (symbolic-simp-rule '(+ (* e_fact (* (+ (* 0.0 c) (* c 0.0)) phi)) (* e_fact (* (* c c) 1.0))))
     (symbolic-simp-rule '(* e_fact (* (+ (* 0.0 c) (* c 0.0)) phi)))
      (symbolic-simp-rule 'e_fact)
      'e_fact
      (symbolic-simp-rule '(* (+ (* 0.0 c) (* c 0.0)) phi))
       (symbolic-simp-rule '(+ (* 0.0 c) (* c 0.0)))
        (symbolic-simp-rule '(* 0.0 c))
        0.0
        (symbolic-simp-rule '(* c 0.0))
        '(* 0.0 c)
       '(+ 0.0 (* 0.0 c))
       (symbolic-simp-rule 'phi)
       'phi
      '(* (+ 0.0 (* 0.0 c)) phi)
     '(* e_fact (* (+ 0.0 (* 0.0 c)) phi))
     (symbolic-simp-rule '(* e_fact (* (* c c) 1.0)))
      (symbolic-simp-rule 'e_fact)
      'e_fact
      (symbolic-simp-rule '(* (* c c) 1.0))
      '(* c (* c 1.0))
     '(* e_fact (* c (* c 1.0)))
    '(+ (* e_fact (* (+ 0.0 (* 0.0 c)) phi)) (* e_fact (* c (* c 1.0))))
   (symbolic-simp '(+ (* e_fact (* (+ 0.0 (* 0.0 c)) phi)) (* e_fact (* c (* c 1.0)))))
    (symbolic-simp-rule '(+ (* e_fact (* (+ 0.0 (* 0.0 c)) phi)) (* e_fact (* c (* c 1.0)))))
     (symbolic-simp-rule '(* e_fact (* (+ 0.0 (* 0.0 c)) phi)))
      (symbolic-simp-rule 'e_fact)
      'e_fact
      (symbolic-simp-rule '(* (+ 0.0 (* 0.0 c)) phi))
       (symbolic-simp-rule '(+ 0.0 (* 0.0 c)))
       '(* 0.0 c)
       (symbolic-simp-rule 'phi)
       'phi
      '(* (* 0.0 c) phi)
     '(* e_fact (* (* 0.0 c) phi))
     (symbolic-simp-rule '(* e_fact (* c (* c 1.0))))
      (symbolic-simp-rule 'e_fact)
      'e_fact
      (symbolic-simp-rule '(* c (* c 1.0)))
       (symbolic-simp-rule 'c)
       'c
       (symbolic-simp-rule '(* c 1.0))
       '(* 1.0 c)
      '(* c (* 1.0 c))
     '(* e_fact (* c (* 1.0 c)))
    '(+ (* e_fact (* (* 0.0 c) phi)) (* e_fact (* c (* 1.0 c))))
   (symbolic-simp '(+ (* e_fact (* (* 0.0 c) phi)) (* e_fact (* c (* 1.0 c)))))
    (symbolic-simp-rule '(+ (* e_fact (* (* 0.0 c) phi)) (* e_fact (* c (* 1.0 c)))))
     (symbolic-simp-rule '(* e_fact (* (* 0.0 c) phi)))
      (symbolic-simp-rule 'e_fact)
      'e_fact
      (symbolic-simp-rule '(* (* 0.0 c) phi))
      '(* 0.0 (* c phi))
     '(* e_fact (* 0.0 (* c phi)))
     (symbolic-simp-rule '(* e_fact (* c (* 1.0 c))))
      (symbolic-simp-rule 'e_fact)
      'e_fact
      (symbolic-simp-rule '(* c (* 1.0 c)))
      '(* 1.0 (* c c))
     '(* e_fact (* 1.0 (* c c)))
    '(+ (* e_fact (* 0.0 (* c phi))) (* e_fact (* 1.0 (* c c))))
   (symbolic-simp '(+ (* e_fact (* 0.0 (* c phi))) (* e_fact (* 1.0 (* c c)))))
    (symbolic-simp-rule '(+ (* e_fact (* 0.0 (* c phi))) (* e_fact (* 1.0 (* c c)))))
     (symbolic-simp-rule '(* e_fact (* 0.0 (* c phi))))
     '(* 0.0 (* e_fact (* c phi)))
     (symbolic-simp-rule '(* e_fact (* 1.0 (* c c))))
     '(* 1.0 (* e_fact (* c c)))
    '(+ (* 0.0 (* e_fact (* c phi))) (* 1.0 (* e_fact (* c c))))
   (symbolic-simp '(+ (* 0.0 (* e_fact (* c phi))) (* 1.0 (* e_fact (* c c)))))
    (symbolic-simp-rule '(+ (* 0.0 (* e_fact (* c phi))) (* 1.0 (* e_fact (* c c)))))
     (symbolic-simp-rule '(* 0.0 (* e_fact (* c phi))))
     0.0
     (symbolic-simp-rule '(* 1.0 (* e_fact (* c c))))
     '(* e_fact (* c c))
    '(+ 0.0 (* e_fact (* c c)))
   (symbolic-simp '(+ 0.0 (* e_fact (* c c))))
    (symbolic-simp-rule '(+ 0.0 (* e_fact (* c c))))
    '(* e_fact (* c c))
   (symbolic-simp '(* e_fact (* c c)))
    (symbolic-simp-rule '(* e_fact (* c c)))
     (symbolic-simp-rule 'e_fact)
     'e_fact
     (symbolic-simp-rule '(* c c))
      (symbolic-simp-rule 'c)
      'c
      (symbolic-simp-rule 'c)
      'c
     '(* c c)
    '(* e_fact (* c c))
   '(* e_fact (* c c))
   (symbolic-diff '(* e_fact Ex) 'Ex)
    (symbolic-diff 'e_fact 'Ex)
    0.0
    (symbolic-diff 'Ex 'Ex)
    1.0
   '(+ (* 0.0 Ex) (* e_fact 1.0))
   (symbolic-simp '(+ (* 0.0 Ex) (* e_fact 1.0)))
    (symbolic-simp-rule '(+ (* 0.0 Ex) (* e_fact 1.0)))
     (symbolic-simp-rule '(* 0.0 Ex))
     0.0
     (symbolic-simp-rule '(* e_fact 1.0))
     '(* 1.0 e_fact)
    '(+ 0.0 (* 1.0 e_fact))
   (symbolic-simp '(+ 0.0 (* 1.0 e_fact)))
    (symbolic-simp-rule '(+ 0.0 (* 1.0 e_fact)))
    '(* 1.0 e_fact)
   (symbolic-simp '(* 1.0 e_fact))
    (symbolic-simp-rule '(* 1.0 e_fact))
    'e_fact
   (symbolic-simp 'e_fact)
    (symbolic-simp-rule 'e_fact)
    'e_fact
   'e_fact
   (symbolic-diff '(* e_fact Ex) 'phi)
    (symbolic-diff 'e_fact 'phi)
    0.0
    (symbolic-diff 'Ex 'phi)
    0.0
   '(+ (* 0.0 Ex) (* e_fact 0.0))
   (symbolic-simp '(+ (* 0.0 Ex) (* e_fact 0.0)))
    (symbolic-simp-rule '(+ (* 0.0 Ex) (* e_fact 0.0)))
     (symbolic-simp-rule '(* 0.0 Ex))
     0.0
     (symbolic-simp-rule '(* e_fact 0.0))
     '(* 0.0 e_fact)
    '(+ 0.0 (* 0.0 e_fact))
   (symbolic-simp '(+ 0.0 (* 0.0 e_fact)))
    (symbolic-simp-rule '(+ 0.0 (* 0.0 e_fact)))
    '(* 0.0 e_fact)
   (symbolic-simp '(* 0.0 e_fact))
    (symbolic-simp-rule '(* 0.0 e_fact))
    0.0
   (symbolic-simp 0.0)
    (symbolic-simp-rule 0.0)
    0.0
   0.0
  '((0.0 (* e_fact (* c c))) (e_fact 0.0))
  (symbolic-roe-matrix '((0.0 (* e_fact (* c c))) (e_fact 0.0)) '(Ex phi))
   (flux-deriv-replace 0.0 'Ex 'ExL)
   0.0
   (flux-deriv-replace 0.0 'phi 'phiL)
   0.0
   (flux-deriv-replace 0.0 'Ex 'ExR)
   0.0
   (flux-deriv-replace 0.0 'phi 'phiR)
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
   (flux-deriv-replace '(* e_fact (* c c)) 'Ex 'ExL)
    (flux-deriv-replace 'e_fact 'Ex 'ExL)
    'e_fact
    (flux-deriv-replace '(* c c) 'Ex 'ExL)
     (flux-deriv-replace 'c 'Ex 'ExL)
     'c
     (flux-deriv-replace 'c 'Ex 'ExL)
     'c
    '(* c c)
   '(* e_fact (* c c))
   (flux-deriv-replace '(* e_fact (* c c)) 'phi 'phiL)
    (flux-deriv-replace 'e_fact 'phi 'phiL)
    'e_fact
    (flux-deriv-replace '(* c c) 'phi 'phiL)
     (flux-deriv-replace 'c 'phi 'phiL)
     'c
     (flux-deriv-replace 'c 'phi 'phiL)
     'c
    '(* c c)
   '(* e_fact (* c c))
   (flux-deriv-replace '(* e_fact (* c c)) 'Ex 'ExR)
    (flux-deriv-replace 'e_fact 'Ex 'ExR)
    'e_fact
    (flux-deriv-replace '(* c c) 'Ex 'ExR)
     (flux-deriv-replace 'c 'Ex 'ExR)
     'c
     (flux-deriv-replace 'c 'Ex 'ExR)
     'c
    '(* c c)
   '(* e_fact (* c c))
   (flux-deriv-replace '(* e_fact (* c c)) 'phi 'phiR)
    (flux-deriv-replace 'e_fact 'phi 'phiR)
    'e_fact
    (flux-deriv-replace '(* c c) 'phi 'phiR)
     (flux-deriv-replace 'c 'phi 'phiR)
     'c
     (flux-deriv-replace 'c 'phi 'phiR)
     'c
    '(* c c)
   '(* e_fact (* c c))
   (symbolic-simp '(+ (* 0.5 (* e_fact (* c c))) (* 0.5 (* e_fact (* c c)))))
    (symbolic-simp-rule '(+ (* 0.5 (* e_fact (* c c))) (* 0.5 (* e_fact (* c c)))))
    '(* (+ 0.5 0.5) (* e_fact (* c c)))
   (symbolic-simp '(* (+ 0.5 0.5) (* e_fact (* c c))))
    (symbolic-simp-rule '(* (+ 0.5 0.5) (* e_fact (* c c))))
     (symbolic-simp-rule '(+ 0.5 0.5))
     1.0
     (symbolic-simp-rule '(* e_fact (* c c)))
      (symbolic-simp-rule 'e_fact)
      'e_fact
      (symbolic-simp-rule '(* c c))
       (symbolic-simp-rule 'c)
       'c
       (symbolic-simp-rule 'c)
       'c
      '(* c c)
     '(* e_fact (* c c))
    '(* 1.0 (* e_fact (* c c)))
   (symbolic-simp '(* 1.0 (* e_fact (* c c))))
    (symbolic-simp-rule '(* 1.0 (* e_fact (* c c))))
    '(* e_fact (* c c))
   (symbolic-simp '(* e_fact (* c c)))
    (symbolic-simp-rule '(* e_fact (* c c)))
     (symbolic-simp-rule 'e_fact)
     'e_fact
     (symbolic-simp-rule '(* c c))
      (symbolic-simp-rule 'c)
      'c
      (symbolic-simp-rule 'c)
      'c
     '(* c c)
    '(* e_fact (* c c))
   '(* e_fact (* c c))
   (flux-deriv-replace 'e_fact 'Ex 'ExL)
   'e_fact
   (flux-deriv-replace 'e_fact 'phi 'phiL)
   'e_fact
   (flux-deriv-replace 'e_fact 'Ex 'ExR)
   'e_fact
   (flux-deriv-replace 'e_fact 'phi 'phiR)
   'e_fact
   (symbolic-simp '(+ (* 0.5 e_fact) (* 0.5 e_fact)))
    (symbolic-simp-rule '(+ (* 0.5 e_fact) (* 0.5 e_fact)))
    '(* (+ 0.5 0.5) e_fact)
   (symbolic-simp '(* (+ 0.5 0.5) e_fact))
    (symbolic-simp-rule '(* (+ 0.5 0.5) e_fact))
     (symbolic-simp-rule '(+ 0.5 0.5))
     1.0
     (symbolic-simp-rule 'e_fact)
     'e_fact
    '(* 1.0 e_fact)
   (symbolic-simp '(* 1.0 e_fact))
    (symbolic-simp-rule '(* 1.0 e_fact))
    'e_fact
   (symbolic-simp 'e_fact)
    (symbolic-simp-rule 'e_fact)
    'e_fact
   'e_fact
   (flux-deriv-replace 0.0 'Ex 'ExL)
   0.0
   (flux-deriv-replace 0.0 'phi 'phiL)
   0.0
   (flux-deriv-replace 0.0 'Ex 'ExR)
   0.0
   (flux-deriv-replace 0.0 'phi 'phiR)
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
  '((0.0 (* e_fact (* c c))) (e_fact 0.0))
  (symbolic-simp '(- ExL ExR))
   (symbolic-simp-rule '(- ExL ExR))
    (symbolic-simp-rule 'ExL)
    'ExL
    (symbolic-simp-rule 'ExR)
    'ExR
   '(- ExL ExR)
  '(- ExL ExR)
  (symbolic-simp '(- phiL phiR))
   (symbolic-simp-rule '(- phiL phiR))
    (symbolic-simp-rule 'phiL)
    'phiL
    (symbolic-simp-rule 'phiR)
    'phiR
   '(- phiL phiR)
  '(- phiL phiR)
  (symbolic-simp '(+ (* 0.0 (- ExL ExR)) (* (* e_fact (* c c)) (- phiL phiR))))
   (symbolic-simp-rule '(+ (* 0.0 (- ExL ExR)) (* (* e_fact (* c c)) (- phiL phiR))))
    (symbolic-simp-rule '(* 0.0 (- ExL ExR)))
    0.0
    (symbolic-simp-rule '(* (* e_fact (* c c)) (- phiL phiR)))
    '(* e_fact (* (* c c) (- phiL phiR)))
   '(+ 0.0 (* e_fact (* (* c c) (- phiL phiR))))
  (symbolic-simp '(+ 0.0 (* e_fact (* (* c c) (- phiL phiR)))))
   (symbolic-simp-rule '(+ 0.0 (* e_fact (* (* c c) (- phiL phiR)))))
   '(* e_fact (* (* c c) (- phiL phiR)))
  (symbolic-simp '(* e_fact (* (* c c) (- phiL phiR))))
   (symbolic-simp-rule '(* e_fact (* (* c c) (- phiL phiR))))
    (symbolic-simp-rule 'e_fact)
    'e_fact
    (symbolic-simp-rule '(* (* c c) (- phiL phiR)))
    '(* c (* c (- phiL phiR)))
   '(* e_fact (* c (* c (- phiL phiR))))
  (symbolic-simp '(* e_fact (* c (* c (- phiL phiR)))))
   (symbolic-simp-rule '(* e_fact (* c (* c (- phiL phiR)))))
    (symbolic-simp-rule 'e_fact)
    'e_fact
    (symbolic-simp-rule '(* c (* c (- phiL phiR))))
     (symbolic-simp-rule 'c)
     'c
     (symbolic-simp-rule '(* c (- phiL phiR)))
      (symbolic-simp-rule 'c)
      'c
      (symbolic-simp-rule '(- phiL phiR))
       (symbolic-simp-rule 'phiL)
       'phiL
       (symbolic-simp-rule 'phiR)
       'phiR
      '(- phiL phiR)
     '(* c (- phiL phiR))
    '(* c (* c (- phiL phiR)))
   '(* e_fact (* c (* c (- phiL phiR))))
  '(* e_fact (* c (* c (- phiL phiR))))
  (symbolic-simp '(+ (* e_fact (- ExL ExR)) (* 0.0 (- phiL phiR))))
   (symbolic-simp-rule '(+ (* e_fact (- ExL ExR)) (* 0.0 (- phiL phiR))))
    (symbolic-simp-rule '(* e_fact (- ExL ExR)))
     (symbolic-simp-rule 'e_fact)
     'e_fact
     (symbolic-simp-rule '(- ExL ExR))
      (symbolic-simp-rule 'ExL)
      'ExL
      (symbolic-simp-rule 'ExR)
      'ExR
     '(- ExL ExR)
    '(* e_fact (- ExL ExR))
    (symbolic-simp-rule '(* 0.0 (- phiL phiR)))
    0.0
   '(+ (* e_fact (- ExL ExR)) 0.0)
  (symbolic-simp '(+ (* e_fact (- ExL ExR)) 0.0))
   (symbolic-simp-rule '(+ (* e_fact (- ExL ExR)) 0.0))
   '(+ 0.0 (* e_fact (- ExL ExR)))
  (symbolic-simp '(+ 0.0 (* e_fact (- ExL ExR))))
   (symbolic-simp-rule '(+ 0.0 (* e_fact (- ExL ExR))))
   '(* e_fact (- ExL ExR))
  (symbolic-simp '(* e_fact (- ExL ExR)))
   (symbolic-simp-rule '(* e_fact (- ExL ExR)))
    (symbolic-simp-rule 'e_fact)
    'e_fact
    (symbolic-simp-rule '(- ExL ExR))
     (symbolic-simp-rule 'ExL)
     'ExL
     (symbolic-simp-rule 'ExR)
     'ExR
    '(- ExL ExR)
   '(* e_fact (- ExL ExR))
  '(* e_fact (- ExL ExR))
  (flux-deriv-replace '(* e_fact (* (* c c) phi)) 'Ex 'ExL)
   (flux-deriv-replace 'e_fact 'Ex 'ExL)
   'e_fact
   (flux-deriv-replace '(* (* c c) phi) 'Ex 'ExL)
    (flux-deriv-replace '(* c c) 'Ex 'ExL)
     (flux-deriv-replace 'c 'Ex 'ExL)
     'c
     (flux-deriv-replace 'c 'Ex 'ExL)
     'c
    '(* c c)
    (flux-deriv-replace 'phi 'Ex 'ExL)
    'phi
   '(* (* c c) phi)
  '(* e_fact (* (* c c) phi))
  (flux-deriv-replace '(* e_fact (* (* c c) phi)) 'phi 'phiL)
   (flux-deriv-replace 'e_fact 'phi 'phiL)
   'e_fact
   (flux-deriv-replace '(* (* c c) phi) 'phi 'phiL)
    (flux-deriv-replace '(* c c) 'phi 'phiL)
     (flux-deriv-replace 'c 'phi 'phiL)
     'c
     (flux-deriv-replace 'c 'phi 'phiL)
     'c
    '(* c c)
    (flux-deriv-replace 'phi 'phi 'phiL)
    'phiL
   '(* (* c c) phiL)
  '(* e_fact (* (* c c) phiL))
  (flux-deriv-replace '(* e_fact (* (* c c) phi)) 'Ex 'ExR)
   (flux-deriv-replace 'e_fact 'Ex 'ExR)
   'e_fact
   (flux-deriv-replace '(* (* c c) phi) 'Ex 'ExR)
    (flux-deriv-replace '(* c c) 'Ex 'ExR)
     (flux-deriv-replace 'c 'Ex 'ExR)
     'c
     (flux-deriv-replace 'c 'Ex 'ExR)
     'c
    '(* c c)
    (flux-deriv-replace 'phi 'Ex 'ExR)
    'phi
   '(* (* c c) phi)
  '(* e_fact (* (* c c) phi))
  (flux-deriv-replace '(* e_fact (* (* c c) phi)) 'phi 'phiR)
   (flux-deriv-replace 'e_fact 'phi 'phiR)
   'e_fact
   (flux-deriv-replace '(* (* c c) phi) 'phi 'phiR)
    (flux-deriv-replace '(* c c) 'phi 'phiR)
     (flux-deriv-replace 'c 'phi 'phiR)
     'c
     (flux-deriv-replace 'c 'phi 'phiR)
     'c
    '(* c c)
    (flux-deriv-replace 'phi 'phi 'phiR)
    'phiR
   '(* (* c c) phiR)
  '(* e_fact (* (* c c) phiR))
  (symbolic-simp '(- (* e_fact (* (* c c) phiL)) (* e_fact (* (* c c) phiR))))
   (symbolic-simp-rule '(- (* e_fact (* (* c c) phiL)) (* e_fact (* (* c c) phiR))))
   '(* e_fact (- (* (* c c) phiL) (* (* c c) phiR)))
  (symbolic-simp '(* e_fact (- (* (* c c) phiL) (* (* c c) phiR))))
   (symbolic-simp-rule '(* e_fact (- (* (* c c) phiL) (* (* c c) phiR))))
    (symbolic-simp-rule 'e_fact)
    'e_fact
    (symbolic-simp-rule '(- (* (* c c) phiL) (* (* c c) phiR)))
    '(* (* c c) (- phiL phiR))
   '(* e_fact (* (* c c) (- phiL phiR)))
  (symbolic-simp '(* e_fact (* (* c c) (- phiL phiR))))
   (symbolic-simp-rule '(* e_fact (* (* c c) (- phiL phiR))))
    (symbolic-simp-rule 'e_fact)
    'e_fact
    (symbolic-simp-rule '(* (* c c) (- phiL phiR)))
    '(* c (* c (- phiL phiR)))
   '(* e_fact (* c (* c (- phiL phiR))))
  (symbolic-simp '(* e_fact (* c (* c (- phiL phiR)))))
   (symbolic-simp-rule '(* e_fact (* c (* c (- phiL phiR)))))
    (symbolic-simp-rule 'e_fact)
    'e_fact
    (symbolic-simp-rule '(* c (* c (- phiL phiR))))
     (symbolic-simp-rule 'c)
     'c
     (symbolic-simp-rule '(* c (- phiL phiR)))
      (symbolic-simp-rule 'c)
      'c
      (symbolic-simp-rule '(- phiL phiR))
       (symbolic-simp-rule 'phiL)
       'phiL
       (symbolic-simp-rule 'phiR)
       'phiR
      '(- phiL phiR)
     '(* c (- phiL phiR))
    '(* c (* c (- phiL phiR)))
   '(* e_fact (* c (* c (- phiL phiR))))
  '(* e_fact (* c (* c (- phiL phiR))))
  (flux-deriv-replace '(* e_fact Ex) 'Ex 'ExL)
   (flux-deriv-replace 'e_fact 'Ex 'ExL)
   'e_fact
   (flux-deriv-replace 'Ex 'Ex 'ExL)
   'ExL
  '(* e_fact ExL)
  (flux-deriv-replace '(* e_fact ExL) 'phi 'phiL)
   (flux-deriv-replace 'e_fact 'phi 'phiL)
   'e_fact
   (flux-deriv-replace 'ExL 'phi 'phiL)
   'ExL
  '(* e_fact ExL)
  (flux-deriv-replace '(* e_fact Ex) 'Ex 'ExR)
   (flux-deriv-replace 'e_fact 'Ex 'ExR)
   'e_fact
   (flux-deriv-replace 'Ex 'Ex 'ExR)
   'ExR
  '(* e_fact ExR)
  (flux-deriv-replace '(* e_fact ExR) 'phi 'phiR)
   (flux-deriv-replace 'e_fact 'phi 'phiR)
   'e_fact
   (flux-deriv-replace 'ExR 'phi 'phiR)
   'ExR
  '(* e_fact ExR)
  (symbolic-simp '(- (* e_fact ExL) (* e_fact ExR)))
   (symbolic-simp-rule '(- (* e_fact ExL) (* e_fact ExR)))
   '(* e_fact (- ExL ExR))
  (symbolic-simp '(* e_fact (- ExL ExR)))
   (symbolic-simp-rule '(* e_fact (- ExL ExR)))
    (symbolic-simp-rule 'e_fact)
    'e_fact
    (symbolic-simp-rule '(- ExL ExR))
     (symbolic-simp-rule 'ExL)
     'ExL
     (symbolic-simp-rule 'ExR)
     'ExR
    '(- ExL ExR)
   '(* e_fact (- ExL ExR))
  '(* e_fact (- ExL ExR))
  (is-real 1.0 '(Ex phi) '((define c 1.0) (define e_fact 1.0) (define b_fact 1.0)))
  #t
  (is-real 1.0 '(Ex phi) '((define c 1.0) (define e_fact 1.0) (define b_fact 1.0)))
  #t
  (is-real 1.0 '(Ex phi) '((define c 1.0) (define e_fact 1.0) (define b_fact 1.0)))
  #t
  (is-real 0.0 '(Ex phi) '((define c 1.0) (define e_fact 1.0) (define b_fact 1.0)))
  #t
  (is-real '(cond ((< x 0.0) 0.5) (else -0.5)) '(Ex phi) '((define c 1.0) (define e_fact 1.0) (define b_fact 1.0)))
   (is-real 0.5 '(Ex phi) '((define c 1.0) (define e_fact 1.0) (define b_fact 1.0)))
   #t
  (is-real -0.5 '(Ex phi) '((define c 1.0) (define e_fact 1.0) (define b_fact 1.0)))
  #t
 #t
