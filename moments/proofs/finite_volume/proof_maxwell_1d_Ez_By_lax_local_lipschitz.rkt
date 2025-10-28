#lang racket

(require "../prover_core.rkt")
(require "../prover_vector.rkt")

 (prove-lax-friedrichs-vector2-1d-local-lipschitz '#hash((cons-exprs . (Ez By)) (flux-exprs . ((* -1.0 (* (* c c) By)) (* -1.0 Ez))) (max-speed-exprs . ((abs c) (abs c))) (name . "maxwell_EzBy") (parameters . ((define c 1.0) (define e_fact 1.0) (define b_fact 1.0)))) #:cfl 0.95 #:init-funcs '(0.0 (cond ((< x 0.0) 0.5) (else -0.5))) #:nx 200 #:t-final 1.0 #:x0 -1.5 #:x1 1.5)
  (symbolic-hessian '(* -1.0 (* (* c c) By)) '(Ez By))
   (symbolic-gradient '(* -1.0 (* (* c c) By)) '(Ez By))
    (symbolic-diff '(* -1.0 (* (* c c) By)) 'Ez)
     (symbolic-diff -1.0 'Ez)
     0.0
     (symbolic-diff '(* (* c c) By) 'Ez)
      (symbolic-diff '(* c c) 'Ez)
       (symbolic-diff 'c 'Ez)
       0.0
       (symbolic-diff 'c 'Ez)
       0.0
      '(+ (* 0.0 c) (* c 0.0))
      (symbolic-diff 'By 'Ez)
      0.0
     '(+ (* (+ (* 0.0 c) (* c 0.0)) By) (* (* c c) 0.0))
    '(+ (* 0.0 (* (* c c) By)) (* -1.0 (+ (* (+ (* 0.0 c) (* c 0.0)) By) (* (* c c) 0.0))))
    (symbolic-simp '(+ (* 0.0 (* (* c c) By)) (* -1.0 (+ (* (+ (* 0.0 c) (* c 0.0)) By) (* (* c c) 0.0)))))
     (symbolic-simp-rule '(+ (* 0.0 (* (* c c) By)) (* -1.0 (+ (* (+ (* 0.0 c) (* c 0.0)) By) (* (* c c) 0.0)))))
      (symbolic-simp-rule '(* 0.0 (* (* c c) By)))
      0.0
      (symbolic-simp-rule '(* -1.0 (+ (* (+ (* 0.0 c) (* c 0.0)) By) (* (* c c) 0.0))))
      '(+ (* -1.0 (* (+ (* 0.0 c) (* c 0.0)) By)) (* -1.0 (* (* c c) 0.0)))
     '(+ 0.0 (+ (* -1.0 (* (+ (* 0.0 c) (* c 0.0)) By)) (* -1.0 (* (* c c) 0.0))))
    (symbolic-simp '(+ 0.0 (+ (* -1.0 (* (+ (* 0.0 c) (* c 0.0)) By)) (* -1.0 (* (* c c) 0.0)))))
     (symbolic-simp-rule '(+ 0.0 (+ (* -1.0 (* (+ (* 0.0 c) (* c 0.0)) By)) (* -1.0 (* (* c c) 0.0)))))
     '(+ (* -1.0 (* (+ (* 0.0 c) (* c 0.0)) By)) (* -1.0 (* (* c c) 0.0)))
    (symbolic-simp '(+ (* -1.0 (* (+ (* 0.0 c) (* c 0.0)) By)) (* -1.0 (* (* c c) 0.0))))
     (symbolic-simp-rule '(+ (* -1.0 (* (+ (* 0.0 c) (* c 0.0)) By)) (* -1.0 (* (* c c) 0.0))))
      (symbolic-simp-rule '(* -1.0 (* (+ (* 0.0 c) (* c 0.0)) By)))
       (symbolic-simp-rule -1.0)
       -1.0
       (symbolic-simp-rule '(* (+ (* 0.0 c) (* c 0.0)) By))
        (symbolic-simp-rule '(+ (* 0.0 c) (* c 0.0)))
         (symbolic-simp-rule '(* 0.0 c))
         0.0
         (symbolic-simp-rule '(* c 0.0))
         '(* 0.0 c)
        '(+ 0.0 (* 0.0 c))
        (symbolic-simp-rule 'By)
        'By
       '(* (+ 0.0 (* 0.0 c)) By)
      '(* -1.0 (* (+ 0.0 (* 0.0 c)) By))
      (symbolic-simp-rule '(* -1.0 (* (* c c) 0.0)))
       (symbolic-simp-rule -1.0)
       -1.0
       (symbolic-simp-rule '(* (* c c) 0.0))
       '(* c (* c 0.0))
      '(* -1.0 (* c (* c 0.0)))
     '(+ (* -1.0 (* (+ 0.0 (* 0.0 c)) By)) (* -1.0 (* c (* c 0.0))))
    (symbolic-simp '(+ (* -1.0 (* (+ 0.0 (* 0.0 c)) By)) (* -1.0 (* c (* c 0.0)))))
     (symbolic-simp-rule '(+ (* -1.0 (* (+ 0.0 (* 0.0 c)) By)) (* -1.0 (* c (* c 0.0)))))
      (symbolic-simp-rule '(* -1.0 (* (+ 0.0 (* 0.0 c)) By)))
       (symbolic-simp-rule -1.0)
       -1.0
       (symbolic-simp-rule '(* (+ 0.0 (* 0.0 c)) By))
        (symbolic-simp-rule '(+ 0.0 (* 0.0 c)))
        '(* 0.0 c)
        (symbolic-simp-rule 'By)
        'By
       '(* (* 0.0 c) By)
      '(* -1.0 (* (* 0.0 c) By))
      (symbolic-simp-rule '(* -1.0 (* c (* c 0.0))))
       (symbolic-simp-rule -1.0)
       -1.0
       (symbolic-simp-rule '(* c (* c 0.0)))
        (symbolic-simp-rule 'c)
        'c
        (symbolic-simp-rule '(* c 0.0))
        '(* 0.0 c)
       '(* c (* 0.0 c))
      '(* -1.0 (* c (* 0.0 c)))
     '(+ (* -1.0 (* (* 0.0 c) By)) (* -1.0 (* c (* 0.0 c))))
    (symbolic-simp '(+ (* -1.0 (* (* 0.0 c) By)) (* -1.0 (* c (* 0.0 c)))))
     (symbolic-simp-rule '(+ (* -1.0 (* (* 0.0 c) By)) (* -1.0 (* c (* 0.0 c)))))
      (symbolic-simp-rule '(* -1.0 (* (* 0.0 c) By)))
       (symbolic-simp-rule -1.0)
       -1.0
       (symbolic-simp-rule '(* (* 0.0 c) By))
       '(* 0.0 (* c By))
      '(* -1.0 (* 0.0 (* c By)))
      (symbolic-simp-rule '(* -1.0 (* c (* 0.0 c))))
       (symbolic-simp-rule -1.0)
       -1.0
       (symbolic-simp-rule '(* c (* 0.0 c)))
       '(* 0.0 (* c c))
      '(* -1.0 (* 0.0 (* c c)))
     '(+ (* -1.0 (* 0.0 (* c By))) (* -1.0 (* 0.0 (* c c))))
    (symbolic-simp '(+ (* -1.0 (* 0.0 (* c By))) (* -1.0 (* 0.0 (* c c)))))
     (symbolic-simp-rule '(+ (* -1.0 (* 0.0 (* c By))) (* -1.0 (* 0.0 (* c c)))))
      (symbolic-simp-rule '(* -1.0 (* 0.0 (* c By))))
      '(* -0.0 (* c By))
      (symbolic-simp-rule '(* -1.0 (* 0.0 (* c c))))
      '(* -0.0 (* c c))
     '(+ (* -0.0 (* c By)) (* -0.0 (* c c)))
    (symbolic-simp '(+ (* -0.0 (* c By)) (* -0.0 (* c c))))
     (symbolic-simp-rule '(+ (* -0.0 (* c By)) (* -0.0 (* c c))))
      (symbolic-simp-rule '(* -0.0 (* c By)))
      0.0
      (symbolic-simp-rule '(* -0.0 (* c c)))
      0.0
     '(+ 0.0 0.0)
    (symbolic-simp '(+ 0.0 0.0))
     (symbolic-simp-rule '(+ 0.0 0.0))
     0.0
    (symbolic-simp 0.0)
     (symbolic-simp-rule 0.0)
     0.0
    0.0
    (symbolic-diff '(* -1.0 (* (* c c) By)) 'By)
     (symbolic-diff -1.0 'By)
     0.0
     (symbolic-diff '(* (* c c) By) 'By)
      (symbolic-diff '(* c c) 'By)
       (symbolic-diff 'c 'By)
       0.0
       (symbolic-diff 'c 'By)
       0.0
      '(+ (* 0.0 c) (* c 0.0))
      (symbolic-diff 'By 'By)
      1.0
     '(+ (* (+ (* 0.0 c) (* c 0.0)) By) (* (* c c) 1.0))
    '(+ (* 0.0 (* (* c c) By)) (* -1.0 (+ (* (+ (* 0.0 c) (* c 0.0)) By) (* (* c c) 1.0))))
    (symbolic-simp '(+ (* 0.0 (* (* c c) By)) (* -1.0 (+ (* (+ (* 0.0 c) (* c 0.0)) By) (* (* c c) 1.0)))))
     (symbolic-simp-rule '(+ (* 0.0 (* (* c c) By)) (* -1.0 (+ (* (+ (* 0.0 c) (* c 0.0)) By) (* (* c c) 1.0)))))
      (symbolic-simp-rule '(* 0.0 (* (* c c) By)))
      0.0
      (symbolic-simp-rule '(* -1.0 (+ (* (+ (* 0.0 c) (* c 0.0)) By) (* (* c c) 1.0))))
      '(+ (* -1.0 (* (+ (* 0.0 c) (* c 0.0)) By)) (* -1.0 (* (* c c) 1.0)))
     '(+ 0.0 (+ (* -1.0 (* (+ (* 0.0 c) (* c 0.0)) By)) (* -1.0 (* (* c c) 1.0))))
    (symbolic-simp '(+ 0.0 (+ (* -1.0 (* (+ (* 0.0 c) (* c 0.0)) By)) (* -1.0 (* (* c c) 1.0)))))
     (symbolic-simp-rule '(+ 0.0 (+ (* -1.0 (* (+ (* 0.0 c) (* c 0.0)) By)) (* -1.0 (* (* c c) 1.0)))))
     '(+ (* -1.0 (* (+ (* 0.0 c) (* c 0.0)) By)) (* -1.0 (* (* c c) 1.0)))
    (symbolic-simp '(+ (* -1.0 (* (+ (* 0.0 c) (* c 0.0)) By)) (* -1.0 (* (* c c) 1.0))))
     (symbolic-simp-rule '(+ (* -1.0 (* (+ (* 0.0 c) (* c 0.0)) By)) (* -1.0 (* (* c c) 1.0))))
      (symbolic-simp-rule '(* -1.0 (* (+ (* 0.0 c) (* c 0.0)) By)))
       (symbolic-simp-rule -1.0)
       -1.0
       (symbolic-simp-rule '(* (+ (* 0.0 c) (* c 0.0)) By))
        (symbolic-simp-rule '(+ (* 0.0 c) (* c 0.0)))
         (symbolic-simp-rule '(* 0.0 c))
         0.0
         (symbolic-simp-rule '(* c 0.0))
         '(* 0.0 c)
        '(+ 0.0 (* 0.0 c))
        (symbolic-simp-rule 'By)
        'By
       '(* (+ 0.0 (* 0.0 c)) By)
      '(* -1.0 (* (+ 0.0 (* 0.0 c)) By))
      (symbolic-simp-rule '(* -1.0 (* (* c c) 1.0)))
       (symbolic-simp-rule -1.0)
       -1.0
       (symbolic-simp-rule '(* (* c c) 1.0))
       '(* c (* c 1.0))
      '(* -1.0 (* c (* c 1.0)))
     '(+ (* -1.0 (* (+ 0.0 (* 0.0 c)) By)) (* -1.0 (* c (* c 1.0))))
    (symbolic-simp '(+ (* -1.0 (* (+ 0.0 (* 0.0 c)) By)) (* -1.0 (* c (* c 1.0)))))
     (symbolic-simp-rule '(+ (* -1.0 (* (+ 0.0 (* 0.0 c)) By)) (* -1.0 (* c (* c 1.0)))))
      (symbolic-simp-rule '(* -1.0 (* (+ 0.0 (* 0.0 c)) By)))
       (symbolic-simp-rule -1.0)
       -1.0
       (symbolic-simp-rule '(* (+ 0.0 (* 0.0 c)) By))
        (symbolic-simp-rule '(+ 0.0 (* 0.0 c)))
        '(* 0.0 c)
        (symbolic-simp-rule 'By)
        'By
       '(* (* 0.0 c) By)
      '(* -1.0 (* (* 0.0 c) By))
      (symbolic-simp-rule '(* -1.0 (* c (* c 1.0))))
       (symbolic-simp-rule -1.0)
       -1.0
       (symbolic-simp-rule '(* c (* c 1.0)))
        (symbolic-simp-rule 'c)
        'c
        (symbolic-simp-rule '(* c 1.0))
        '(* 1.0 c)
       '(* c (* 1.0 c))
      '(* -1.0 (* c (* 1.0 c)))
     '(+ (* -1.0 (* (* 0.0 c) By)) (* -1.0 (* c (* 1.0 c))))
    (symbolic-simp '(+ (* -1.0 (* (* 0.0 c) By)) (* -1.0 (* c (* 1.0 c)))))
     (symbolic-simp-rule '(+ (* -1.0 (* (* 0.0 c) By)) (* -1.0 (* c (* 1.0 c)))))
      (symbolic-simp-rule '(* -1.0 (* (* 0.0 c) By)))
       (symbolic-simp-rule -1.0)
       -1.0
       (symbolic-simp-rule '(* (* 0.0 c) By))
       '(* 0.0 (* c By))
      '(* -1.0 (* 0.0 (* c By)))
      (symbolic-simp-rule '(* -1.0 (* c (* 1.0 c))))
       (symbolic-simp-rule -1.0)
       -1.0
       (symbolic-simp-rule '(* c (* 1.0 c)))
       '(* 1.0 (* c c))
      '(* -1.0 (* 1.0 (* c c)))
     '(+ (* -1.0 (* 0.0 (* c By))) (* -1.0 (* 1.0 (* c c))))
    (symbolic-simp '(+ (* -1.0 (* 0.0 (* c By))) (* -1.0 (* 1.0 (* c c)))))
     (symbolic-simp-rule '(+ (* -1.0 (* 0.0 (* c By))) (* -1.0 (* 1.0 (* c c)))))
      (symbolic-simp-rule '(* -1.0 (* 0.0 (* c By))))
      '(* -0.0 (* c By))
      (symbolic-simp-rule '(* -1.0 (* 1.0 (* c c))))
      '(* -1.0 (* c c))
     '(+ (* -0.0 (* c By)) (* -1.0 (* c c)))
    (symbolic-simp '(+ (* -0.0 (* c By)) (* -1.0 (* c c))))
     (symbolic-simp-rule '(+ (* -0.0 (* c By)) (* -1.0 (* c c))))
      (symbolic-simp-rule '(* -0.0 (* c By)))
      0.0
      (symbolic-simp-rule '(* -1.0 (* c c)))
       (symbolic-simp-rule -1.0)
       -1.0
       (symbolic-simp-rule '(* c c))
        (symbolic-simp-rule 'c)
        'c
        (symbolic-simp-rule 'c)
        'c
       '(* c c)
      '(* -1.0 (* c c))
     '(+ 0.0 (* -1.0 (* c c)))
    (symbolic-simp '(+ 0.0 (* -1.0 (* c c))))
     (symbolic-simp-rule '(+ 0.0 (* -1.0 (* c c))))
     '(* -1.0 (* c c))
    (symbolic-simp '(* -1.0 (* c c)))
     (symbolic-simp-rule '(* -1.0 (* c c)))
      (symbolic-simp-rule -1.0)
      -1.0
      (symbolic-simp-rule '(* c c))
       (symbolic-simp-rule 'c)
       'c
       (symbolic-simp-rule 'c)
       'c
      '(* c c)
     '(* -1.0 (* c c))
    '(* -1.0 (* c c))
   '(0.0 (* -1.0 (* c c)))
  (symbolic-jacobian '(0.0 (* -1.0 (* c c))) '(Ez By))
   (symbolic-diff 0.0 'Ez)
   0.0
   (symbolic-simp 0.0)
    (symbolic-simp-rule 0.0)
    0.0
   0.0
   (symbolic-diff 0.0 'By)
   0.0
   (symbolic-simp 0.0)
    (symbolic-simp-rule 0.0)
    0.0
   0.0
   (symbolic-diff '(* -1.0 (* c c)) 'Ez)
    (symbolic-diff -1.0 'Ez)
    0.0
    (symbolic-diff '(* c c) 'Ez)
     (symbolic-diff 'c 'Ez)
     0.0
     (symbolic-diff 'c 'Ez)
     0.0
    '(+ (* 0.0 c) (* c 0.0))
   '(+ (* 0.0 (* c c)) (* -1.0 (+ (* 0.0 c) (* c 0.0))))
   (symbolic-simp '(+ (* 0.0 (* c c)) (* -1.0 (+ (* 0.0 c) (* c 0.0)))))
    (symbolic-simp-rule '(+ (* 0.0 (* c c)) (* -1.0 (+ (* 0.0 c) (* c 0.0)))))
     (symbolic-simp-rule '(* 0.0 (* c c)))
     0.0
     (symbolic-simp-rule '(* -1.0 (+ (* 0.0 c) (* c 0.0))))
     '(+ (* -1.0 (* 0.0 c)) (* -1.0 (* c 0.0)))
    '(+ 0.0 (+ (* -1.0 (* 0.0 c)) (* -1.0 (* c 0.0))))
   (symbolic-simp '(+ 0.0 (+ (* -1.0 (* 0.0 c)) (* -1.0 (* c 0.0)))))
    (symbolic-simp-rule '(+ 0.0 (+ (* -1.0 (* 0.0 c)) (* -1.0 (* c 0.0)))))
    '(+ (* -1.0 (* 0.0 c)) (* -1.0 (* c 0.0)))
   (symbolic-simp '(+ (* -1.0 (* 0.0 c)) (* -1.0 (* c 0.0))))
    (symbolic-simp-rule '(+ (* -1.0 (* 0.0 c)) (* -1.0 (* c 0.0))))
     (symbolic-simp-rule '(* -1.0 (* 0.0 c)))
     '(* -0.0 c)
     (symbolic-simp-rule '(* -1.0 (* c 0.0)))
      (symbolic-simp-rule -1.0)
      -1.0
      (symbolic-simp-rule '(* c 0.0))
      '(* 0.0 c)
     '(* -1.0 (* 0.0 c))
    '(+ (* -0.0 c) (* -1.0 (* 0.0 c)))
   (symbolic-simp '(+ (* -0.0 c) (* -1.0 (* 0.0 c))))
    (symbolic-simp-rule '(+ (* -0.0 c) (* -1.0 (* 0.0 c))))
     (symbolic-simp-rule '(* -0.0 c))
     0.0
     (symbolic-simp-rule '(* -1.0 (* 0.0 c)))
     '(* -0.0 c)
    '(+ 0.0 (* -0.0 c))
   (symbolic-simp '(+ 0.0 (* -0.0 c)))
    (symbolic-simp-rule '(+ 0.0 (* -0.0 c)))
    '(* -0.0 c)
   (symbolic-simp '(* -0.0 c))
    (symbolic-simp-rule '(* -0.0 c))
    0.0
   (symbolic-simp 0.0)
    (symbolic-simp-rule 0.0)
    0.0
   0.0
   (symbolic-diff '(* -1.0 (* c c)) 'By)
    (symbolic-diff -1.0 'By)
    0.0
    (symbolic-diff '(* c c) 'By)
     (symbolic-diff 'c 'By)
     0.0
     (symbolic-diff 'c 'By)
     0.0
    '(+ (* 0.0 c) (* c 0.0))
   '(+ (* 0.0 (* c c)) (* -1.0 (+ (* 0.0 c) (* c 0.0))))
   (symbolic-simp '(+ (* 0.0 (* c c)) (* -1.0 (+ (* 0.0 c) (* c 0.0)))))
    (symbolic-simp-rule '(+ (* 0.0 (* c c)) (* -1.0 (+ (* 0.0 c) (* c 0.0)))))
     (symbolic-simp-rule '(* 0.0 (* c c)))
     0.0
     (symbolic-simp-rule '(* -1.0 (+ (* 0.0 c) (* c 0.0))))
     '(+ (* -1.0 (* 0.0 c)) (* -1.0 (* c 0.0)))
    '(+ 0.0 (+ (* -1.0 (* 0.0 c)) (* -1.0 (* c 0.0))))
   (symbolic-simp '(+ 0.0 (+ (* -1.0 (* 0.0 c)) (* -1.0 (* c 0.0)))))
    (symbolic-simp-rule '(+ 0.0 (+ (* -1.0 (* 0.0 c)) (* -1.0 (* c 0.0)))))
    '(+ (* -1.0 (* 0.0 c)) (* -1.0 (* c 0.0)))
   (symbolic-simp '(+ (* -1.0 (* 0.0 c)) (* -1.0 (* c 0.0))))
    (symbolic-simp-rule '(+ (* -1.0 (* 0.0 c)) (* -1.0 (* c 0.0))))
     (symbolic-simp-rule '(* -1.0 (* 0.0 c)))
     '(* -0.0 c)
     (symbolic-simp-rule '(* -1.0 (* c 0.0)))
      (symbolic-simp-rule -1.0)
      -1.0
      (symbolic-simp-rule '(* c 0.0))
      '(* 0.0 c)
     '(* -1.0 (* 0.0 c))
    '(+ (* -0.0 c) (* -1.0 (* 0.0 c)))
   (symbolic-simp '(+ (* -0.0 c) (* -1.0 (* 0.0 c))))
    (symbolic-simp-rule '(+ (* -0.0 c) (* -1.0 (* 0.0 c))))
     (symbolic-simp-rule '(* -0.0 c))
     0.0
     (symbolic-simp-rule '(* -1.0 (* 0.0 c)))
     '(* -0.0 c)
    '(+ 0.0 (* -0.0 c))
   (symbolic-simp '(+ 0.0 (* -0.0 c)))
    (symbolic-simp-rule '(+ 0.0 (* -0.0 c)))
    '(* -0.0 c)
   (symbolic-simp '(* -0.0 c))
    (symbolic-simp-rule '(* -0.0 c))
    0.0
   (symbolic-simp 0.0)
    (symbolic-simp-rule 0.0)
    0.0
   0.0
  '((0.0 0.0) (0.0 0.0))
  (symbolic-hessian '(* -1.0 Ez) '(Ez By))
   (symbolic-gradient '(* -1.0 Ez) '(Ez By))
    (symbolic-diff '(* -1.0 Ez) 'Ez)
     (symbolic-diff -1.0 'Ez)
     0.0
     (symbolic-diff 'Ez 'Ez)
     1.0
    '(+ (* 0.0 Ez) (* -1.0 1.0))
    (symbolic-simp '(+ (* 0.0 Ez) (* -1.0 1.0)))
     (symbolic-simp-rule '(+ (* 0.0 Ez) (* -1.0 1.0)))
      (symbolic-simp-rule '(* 0.0 Ez))
      0.0
      (symbolic-simp-rule '(* -1.0 1.0))
      -1.0
     '(+ 0.0 -1.0)
    (symbolic-simp '(+ 0.0 -1.0))
     (symbolic-simp-rule '(+ 0.0 -1.0))
     -1.0
    (symbolic-simp -1.0)
     (symbolic-simp-rule -1.0)
     -1.0
    -1.0
    (symbolic-diff '(* -1.0 Ez) 'By)
     (symbolic-diff -1.0 'By)
     0.0
     (symbolic-diff 'Ez 'By)
     0.0
    '(+ (* 0.0 Ez) (* -1.0 0.0))
    (symbolic-simp '(+ (* 0.0 Ez) (* -1.0 0.0)))
     (symbolic-simp-rule '(+ (* 0.0 Ez) (* -1.0 0.0)))
      (symbolic-simp-rule '(* 0.0 Ez))
      0.0
      (symbolic-simp-rule '(* -1.0 0.0))
      -0.0
     '(+ 0.0 -0.0)
    (symbolic-simp '(+ 0.0 -0.0))
     (symbolic-simp-rule '(+ 0.0 -0.0))
     -0.0
    (symbolic-simp -0.0)
     (symbolic-simp-rule -0.0)
     -0.0
    -0.0
   '(-1.0 -0.0)
  (symbolic-jacobian '(-1.0 -0.0) '(Ez By))
   (symbolic-diff -1.0 'Ez)
   0.0
   (symbolic-simp 0.0)
    (symbolic-simp-rule 0.0)
    0.0
   0.0
   (symbolic-diff -1.0 'By)
   0.0
   (symbolic-simp 0.0)
    (symbolic-simp-rule 0.0)
    0.0
   0.0
   (symbolic-diff -0.0 'Ez)
   0.0
   (symbolic-simp 0.0)
    (symbolic-simp-rule 0.0)
    0.0
   0.0
   (symbolic-diff -0.0 'By)
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
  (is-real 1.0 '(Ez By) '((define c 1.0) (define e_fact 1.0) (define b_fact 1.0)))
  #t
  (is-real 1.0 '(Ez By) '((define c 1.0) (define e_fact 1.0) (define b_fact 1.0)))
  #t
  (is-real 1.0 '(Ez By) '((define c 1.0) (define e_fact 1.0) (define b_fact 1.0)))
  #t
  (is-real 0.0 '(Ez By) '((define c 1.0) (define e_fact 1.0) (define b_fact 1.0)))
  #t
  (is-real '(cond ((< x 0.0) 0.5) (else -0.5)) '(Ez By) '((define c 1.0) (define e_fact 1.0) (define b_fact 1.0)))
   (is-real 0.5 '(Ez By) '((define c 1.0) (define e_fact 1.0) (define b_fact 1.0)))
   #t
  (is-real -0.5 '(Ez By) '((define c 1.0) (define e_fact 1.0) (define b_fact 1.0)))
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
