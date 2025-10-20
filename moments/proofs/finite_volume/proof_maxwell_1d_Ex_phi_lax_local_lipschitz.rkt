#lang racket

(require "../prover_core.rkt")
(require "../prover_vector.rkt")

 (prove-lax-friedrichs-vector2-1d-local-lipschitz '#hash((cons-exprs . (Ex phi)) (flux-exprs . ((* e_fact (* (* c c) phi)) (* e_fact Ex))) (max-speed-exprs . ((abs (* c e_fact)) (abs (* c e_fact)))) (name . "maxwell_Exphi") (parameters . ((define c 1.0) (define e_fact 1.0) (define b_fact 1.0)))) #:cfl 0.95 #:init-funcs '(0.0 (cond ((< x 0.0) 0.5) (else -0.5))) #:nx 200 #:t-final 1.0 #:x0 -1.5 #:x1 1.5)
  (symbolic-hessian '(* e_fact (* (* c c) phi)) '(Ex phi))
   (symbolic-gradient '(* e_fact (* (* c c) phi)) '(Ex phi))
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
   '(0.0 (* e_fact (* c c)))
  (symbolic-jacobian '(0.0 (* e_fact (* c c))) '(Ex phi))
   (symbolic-diff 0.0 'Ex)
   0.0
   (symbolic-simp 0.0)
    (symbolic-simp-rule 0.0)
    0.0
   0.0
   (symbolic-diff 0.0 'phi)
   0.0
   (symbolic-simp 0.0)
    (symbolic-simp-rule 0.0)
    0.0
   0.0
   (symbolic-diff '(* e_fact (* c c)) 'Ex)
    (symbolic-diff 'e_fact 'Ex)
    0.0
    (symbolic-diff '(* c c) 'Ex)
     (symbolic-diff 'c 'Ex)
     0.0
     (symbolic-diff 'c 'Ex)
     0.0
    '(+ (* 0.0 c) (* c 0.0))
   '(+ (* 0.0 (* c c)) (* e_fact (+ (* 0.0 c) (* c 0.0))))
   (symbolic-simp '(+ (* 0.0 (* c c)) (* e_fact (+ (* 0.0 c) (* c 0.0)))))
    (symbolic-simp-rule '(+ (* 0.0 (* c c)) (* e_fact (+ (* 0.0 c) (* c 0.0)))))
     (symbolic-simp-rule '(* 0.0 (* c c)))
     0.0
     (symbolic-simp-rule '(* e_fact (+ (* 0.0 c) (* c 0.0))))
     '(+ (* e_fact (* 0.0 c)) (* e_fact (* c 0.0)))
    '(+ 0.0 (+ (* e_fact (* 0.0 c)) (* e_fact (* c 0.0))))
   (symbolic-simp '(+ 0.0 (+ (* e_fact (* 0.0 c)) (* e_fact (* c 0.0)))))
    (symbolic-simp-rule '(+ 0.0 (+ (* e_fact (* 0.0 c)) (* e_fact (* c 0.0)))))
    '(+ (* e_fact (* 0.0 c)) (* e_fact (* c 0.0)))
   (symbolic-simp '(+ (* e_fact (* 0.0 c)) (* e_fact (* c 0.0))))
    (symbolic-simp-rule '(+ (* e_fact (* 0.0 c)) (* e_fact (* c 0.0))))
     (symbolic-simp-rule '(* e_fact (* 0.0 c)))
     '(* 0.0 (* e_fact c))
     (symbolic-simp-rule '(* e_fact (* c 0.0)))
      (symbolic-simp-rule 'e_fact)
      'e_fact
      (symbolic-simp-rule '(* c 0.0))
      '(* 0.0 c)
     '(* e_fact (* 0.0 c))
    '(+ (* 0.0 (* e_fact c)) (* e_fact (* 0.0 c)))
   (symbolic-simp '(+ (* 0.0 (* e_fact c)) (* e_fact (* 0.0 c))))
    (symbolic-simp-rule '(+ (* 0.0 (* e_fact c)) (* e_fact (* 0.0 c))))
     (symbolic-simp-rule '(* 0.0 (* e_fact c)))
     0.0
     (symbolic-simp-rule '(* e_fact (* 0.0 c)))
     '(* 0.0 (* e_fact c))
    '(+ 0.0 (* 0.0 (* e_fact c)))
   (symbolic-simp '(+ 0.0 (* 0.0 (* e_fact c))))
    (symbolic-simp-rule '(+ 0.0 (* 0.0 (* e_fact c))))
    '(* 0.0 (* e_fact c))
   (symbolic-simp '(* 0.0 (* e_fact c)))
    (symbolic-simp-rule '(* 0.0 (* e_fact c)))
    0.0
   (symbolic-simp 0.0)
    (symbolic-simp-rule 0.0)
    0.0
   0.0
   (symbolic-diff '(* e_fact (* c c)) 'phi)
    (symbolic-diff 'e_fact 'phi)
    0.0
    (symbolic-diff '(* c c) 'phi)
     (symbolic-diff 'c 'phi)
     0.0
     (symbolic-diff 'c 'phi)
     0.0
    '(+ (* 0.0 c) (* c 0.0))
   '(+ (* 0.0 (* c c)) (* e_fact (+ (* 0.0 c) (* c 0.0))))
   (symbolic-simp '(+ (* 0.0 (* c c)) (* e_fact (+ (* 0.0 c) (* c 0.0)))))
    (symbolic-simp-rule '(+ (* 0.0 (* c c)) (* e_fact (+ (* 0.0 c) (* c 0.0)))))
     (symbolic-simp-rule '(* 0.0 (* c c)))
     0.0
     (symbolic-simp-rule '(* e_fact (+ (* 0.0 c) (* c 0.0))))
     '(+ (* e_fact (* 0.0 c)) (* e_fact (* c 0.0)))
    '(+ 0.0 (+ (* e_fact (* 0.0 c)) (* e_fact (* c 0.0))))
   (symbolic-simp '(+ 0.0 (+ (* e_fact (* 0.0 c)) (* e_fact (* c 0.0)))))
    (symbolic-simp-rule '(+ 0.0 (+ (* e_fact (* 0.0 c)) (* e_fact (* c 0.0)))))
    '(+ (* e_fact (* 0.0 c)) (* e_fact (* c 0.0)))
   (symbolic-simp '(+ (* e_fact (* 0.0 c)) (* e_fact (* c 0.0))))
    (symbolic-simp-rule '(+ (* e_fact (* 0.0 c)) (* e_fact (* c 0.0))))
     (symbolic-simp-rule '(* e_fact (* 0.0 c)))
     '(* 0.0 (* e_fact c))
     (symbolic-simp-rule '(* e_fact (* c 0.0)))
      (symbolic-simp-rule 'e_fact)
      'e_fact
      (symbolic-simp-rule '(* c 0.0))
      '(* 0.0 c)
     '(* e_fact (* 0.0 c))
    '(+ (* 0.0 (* e_fact c)) (* e_fact (* 0.0 c)))
   (symbolic-simp '(+ (* 0.0 (* e_fact c)) (* e_fact (* 0.0 c))))
    (symbolic-simp-rule '(+ (* 0.0 (* e_fact c)) (* e_fact (* 0.0 c))))
     (symbolic-simp-rule '(* 0.0 (* e_fact c)))
     0.0
     (symbolic-simp-rule '(* e_fact (* 0.0 c)))
     '(* 0.0 (* e_fact c))
    '(+ 0.0 (* 0.0 (* e_fact c)))
   (symbolic-simp '(+ 0.0 (* 0.0 (* e_fact c))))
    (symbolic-simp-rule '(+ 0.0 (* 0.0 (* e_fact c))))
    '(* 0.0 (* e_fact c))
   (symbolic-simp '(* 0.0 (* e_fact c)))
    (symbolic-simp-rule '(* 0.0 (* e_fact c)))
    0.0
   (symbolic-simp 0.0)
    (symbolic-simp-rule 0.0)
    0.0
   0.0
  '((0.0 0.0) (0.0 0.0))
  (symbolic-hessian '(* e_fact Ex) '(Ex phi))
   (symbolic-gradient '(* e_fact Ex) '(Ex phi))
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
   '(e_fact 0.0)
  (symbolic-jacobian '(e_fact 0.0) '(Ex phi))
   (symbolic-diff 'e_fact 'Ex)
   0.0
   (symbolic-simp 0.0)
    (symbolic-simp-rule 0.0)
    0.0
   0.0
   (symbolic-diff 'e_fact 'phi)
   0.0
   (symbolic-simp 0.0)
    (symbolic-simp-rule 0.0)
    0.0
   0.0
   (symbolic-diff 0.0 'Ex)
   0.0
   (symbolic-simp 0.0)
    (symbolic-simp-rule 0.0)
    0.0
   0.0
   (symbolic-diff 0.0 'phi)
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
  (is-non-negative 0.0 '((define c 1.0) (define e_fact 1.0) (define b_fact 1.0)))
  #t
  (is-non-negative 0.0 '((define c 1.0) (define e_fact 1.0) (define b_fact 1.0)))
  #t
  (is-non-negative 0.0 '((define c 1.0) (define e_fact 1.0) (define b_fact 1.0)))
  #t
  (is-non-negative 0.0 '((define c 1.0) (define e_fact 1.0) (define b_fact 1.0)))
  #t
 #t
