#lang racket

(require "../prover_core.rkt")
(require "../prover_vector.rkt")

 (prove-roe-vector2-1d-flux-conservation '#hash((cons-exprs . (Ez By)) (flux-exprs . ((* -1.0 (* (* c c) By)) (* -1.0 Ez))) (max-speed-exprs . ((abs c) (abs c))) (name . "maxwell_EzBy") (parameters . ((define c 1.0) (define e_fact 1.0) (define b_fact 1.0)))) #:cfl 0.95 #:init-funcs '(0.0 (cond ((< x 0.0) 0.5) (else -0.5))) #:nx 200 #:t-final 1.0 #:x0 -1.5 #:x1 1.5)
  (symbolic-jacobian '((* -1.0 (* (* c c) By)) (* -1.0 Ez)) '(Ez By))
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
  '((0.0 (* -1.0 (* c c))) (-1.0 -0.0))
  (symbolic-roe-matrix '((0.0 (* -1.0 (* c c))) (-1.0 -0.0)) '(Ez By))
   (flux-deriv-replace 0.0 'Ez 'EzL)
   0.0
   (flux-deriv-replace 0.0 'By 'ByL)
   0.0
   (flux-deriv-replace 0.0 'Ez 'EzR)
   0.0
   (flux-deriv-replace 0.0 'By 'ByR)
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
   (flux-deriv-replace '(* -1.0 (* c c)) 'Ez 'EzL)
    (flux-deriv-replace -1.0 'Ez 'EzL)
    -1.0
    (flux-deriv-replace '(* c c) 'Ez 'EzL)
     (flux-deriv-replace 'c 'Ez 'EzL)
     'c
     (flux-deriv-replace 'c 'Ez 'EzL)
     'c
    '(* c c)
   '(* -1.0 (* c c))
   (flux-deriv-replace '(* -1.0 (* c c)) 'By 'ByL)
    (flux-deriv-replace -1.0 'By 'ByL)
    -1.0
    (flux-deriv-replace '(* c c) 'By 'ByL)
     (flux-deriv-replace 'c 'By 'ByL)
     'c
     (flux-deriv-replace 'c 'By 'ByL)
     'c
    '(* c c)
   '(* -1.0 (* c c))
   (flux-deriv-replace '(* -1.0 (* c c)) 'Ez 'EzR)
    (flux-deriv-replace -1.0 'Ez 'EzR)
    -1.0
    (flux-deriv-replace '(* c c) 'Ez 'EzR)
     (flux-deriv-replace 'c 'Ez 'EzR)
     'c
     (flux-deriv-replace 'c 'Ez 'EzR)
     'c
    '(* c c)
   '(* -1.0 (* c c))
   (flux-deriv-replace '(* -1.0 (* c c)) 'By 'ByR)
    (flux-deriv-replace -1.0 'By 'ByR)
    -1.0
    (flux-deriv-replace '(* c c) 'By 'ByR)
     (flux-deriv-replace 'c 'By 'ByR)
     'c
     (flux-deriv-replace 'c 'By 'ByR)
     'c
    '(* c c)
   '(* -1.0 (* c c))
   (symbolic-simp '(+ (* 0.5 (* -1.0 (* c c))) (* 0.5 (* -1.0 (* c c)))))
    (symbolic-simp-rule '(+ (* 0.5 (* -1.0 (* c c))) (* 0.5 (* -1.0 (* c c)))))
    '(* (+ 0.5 0.5) (* -1.0 (* c c)))
   (symbolic-simp '(* (+ 0.5 0.5) (* -1.0 (* c c))))
    (symbolic-simp-rule '(* (+ 0.5 0.5) (* -1.0 (* c c))))
    '(* -1.0 (* (+ 0.5 0.5) (* c c)))
   (symbolic-simp '(* -1.0 (* (+ 0.5 0.5) (* c c))))
    (symbolic-simp-rule '(* -1.0 (* (+ 0.5 0.5) (* c c))))
     (symbolic-simp-rule -1.0)
     -1.0
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
    '(* -1.0 (* 1.0 (* c c)))
   (symbolic-simp '(* -1.0 (* 1.0 (* c c))))
    (symbolic-simp-rule '(* -1.0 (* 1.0 (* c c))))
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
   (flux-deriv-replace -1.0 'Ez 'EzL)
   -1.0
   (flux-deriv-replace -1.0 'By 'ByL)
   -1.0
   (flux-deriv-replace -1.0 'Ez 'EzR)
   -1.0
   (flux-deriv-replace -1.0 'By 'ByR)
   -1.0
   (symbolic-simp '(+ (* 0.5 -1.0) (* 0.5 -1.0)))
    (symbolic-simp-rule '(+ (* 0.5 -1.0) (* 0.5 -1.0)))
    '(* (+ 0.5 0.5) -1.0)
   (symbolic-simp '(* (+ 0.5 0.5) -1.0))
    (symbolic-simp-rule '(* (+ 0.5 0.5) -1.0))
    '(* -1.0 (+ 0.5 0.5))
   (symbolic-simp '(* -1.0 (+ 0.5 0.5)))
    (symbolic-simp-rule '(* -1.0 (+ 0.5 0.5)))
    -1.0
   (symbolic-simp -1.0)
    (symbolic-simp-rule -1.0)
    -1.0
   -1.0
   (flux-deriv-replace -0.0 'Ez 'EzL)
   -0.0
   (flux-deriv-replace -0.0 'By 'ByL)
   -0.0
   (flux-deriv-replace -0.0 'Ez 'EzR)
   -0.0
   (flux-deriv-replace -0.0 'By 'ByR)
   -0.0
   (symbolic-simp '(+ (* 0.5 -0.0) (* 0.5 -0.0)))
    (symbolic-simp-rule '(+ (* 0.5 -0.0) (* 0.5 -0.0)))
    '(* (+ 0.5 0.5) -0.0)
   (symbolic-simp '(* (+ 0.5 0.5) -0.0))
    (symbolic-simp-rule '(* (+ 0.5 0.5) -0.0))
    '(* -0.0 (+ 0.5 0.5))
   (symbolic-simp '(* -0.0 (+ 0.5 0.5)))
    (symbolic-simp-rule '(* -0.0 (+ 0.5 0.5)))
    0.0
   (symbolic-simp 0.0)
    (symbolic-simp-rule 0.0)
    0.0
   0.0
  '((0.0 (* -1.0 (* c c))) (-1.0 0.0))
  (symbolic-simp '(- EzL EzR))
   (symbolic-simp-rule '(- EzL EzR))
    (symbolic-simp-rule 'EzL)
    'EzL
    (symbolic-simp-rule 'EzR)
    'EzR
   '(- EzL EzR)
  '(- EzL EzR)
  (symbolic-simp '(- ByL ByR))
   (symbolic-simp-rule '(- ByL ByR))
    (symbolic-simp-rule 'ByL)
    'ByL
    (symbolic-simp-rule 'ByR)
    'ByR
   '(- ByL ByR)
  '(- ByL ByR)
  (symbolic-simp '(+ (* 0.0 (- EzL EzR)) (* (* -1.0 (* c c)) (- ByL ByR))))
   (symbolic-simp-rule '(+ (* 0.0 (- EzL EzR)) (* (* -1.0 (* c c)) (- ByL ByR))))
    (symbolic-simp-rule '(* 0.0 (- EzL EzR)))
    0.0
    (symbolic-simp-rule '(* (* -1.0 (* c c)) (- ByL ByR)))
    '(* -1.0 (* (* c c) (- ByL ByR)))
   '(+ 0.0 (* -1.0 (* (* c c) (- ByL ByR))))
  (symbolic-simp '(+ 0.0 (* -1.0 (* (* c c) (- ByL ByR)))))
   (symbolic-simp-rule '(+ 0.0 (* -1.0 (* (* c c) (- ByL ByR)))))
   '(* -1.0 (* (* c c) (- ByL ByR)))
  (symbolic-simp '(* -1.0 (* (* c c) (- ByL ByR))))
   (symbolic-simp-rule '(* -1.0 (* (* c c) (- ByL ByR))))
    (symbolic-simp-rule -1.0)
    -1.0
    (symbolic-simp-rule '(* (* c c) (- ByL ByR)))
    '(* c (* c (- ByL ByR)))
   '(* -1.0 (* c (* c (- ByL ByR))))
  (symbolic-simp '(* -1.0 (* c (* c (- ByL ByR)))))
   (symbolic-simp-rule '(* -1.0 (* c (* c (- ByL ByR)))))
    (symbolic-simp-rule -1.0)
    -1.0
    (symbolic-simp-rule '(* c (* c (- ByL ByR))))
     (symbolic-simp-rule 'c)
     'c
     (symbolic-simp-rule '(* c (- ByL ByR)))
      (symbolic-simp-rule 'c)
      'c
      (symbolic-simp-rule '(- ByL ByR))
       (symbolic-simp-rule 'ByL)
       'ByL
       (symbolic-simp-rule 'ByR)
       'ByR
      '(- ByL ByR)
     '(* c (- ByL ByR))
    '(* c (* c (- ByL ByR)))
   '(* -1.0 (* c (* c (- ByL ByR))))
  '(* -1.0 (* c (* c (- ByL ByR))))
  (symbolic-simp '(+ (* -1.0 (- EzL EzR)) (* 0.0 (- ByL ByR))))
   (symbolic-simp-rule '(+ (* -1.0 (- EzL EzR)) (* 0.0 (- ByL ByR))))
    (symbolic-simp-rule '(* -1.0 (- EzL EzR)))
     (symbolic-simp-rule -1.0)
     -1.0
     (symbolic-simp-rule '(- EzL EzR))
      (symbolic-simp-rule 'EzL)
      'EzL
      (symbolic-simp-rule 'EzR)
      'EzR
     '(- EzL EzR)
    '(* -1.0 (- EzL EzR))
    (symbolic-simp-rule '(* 0.0 (- ByL ByR)))
    0.0
   '(+ (* -1.0 (- EzL EzR)) 0.0)
  (symbolic-simp '(+ (* -1.0 (- EzL EzR)) 0.0))
   (symbolic-simp-rule '(+ (* -1.0 (- EzL EzR)) 0.0))
   '(+ 0.0 (* -1.0 (- EzL EzR)))
  (symbolic-simp '(+ 0.0 (* -1.0 (- EzL EzR))))
   (symbolic-simp-rule '(+ 0.0 (* -1.0 (- EzL EzR))))
   '(* -1.0 (- EzL EzR))
  (symbolic-simp '(* -1.0 (- EzL EzR)))
   (symbolic-simp-rule '(* -1.0 (- EzL EzR)))
    (symbolic-simp-rule -1.0)
    -1.0
    (symbolic-simp-rule '(- EzL EzR))
     (symbolic-simp-rule 'EzL)
     'EzL
     (symbolic-simp-rule 'EzR)
     'EzR
    '(- EzL EzR)
   '(* -1.0 (- EzL EzR))
  '(* -1.0 (- EzL EzR))
  (flux-deriv-replace '(* -1.0 (* (* c c) By)) 'Ez 'EzL)
   (flux-deriv-replace -1.0 'Ez 'EzL)
   -1.0
   (flux-deriv-replace '(* (* c c) By) 'Ez 'EzL)
    (flux-deriv-replace '(* c c) 'Ez 'EzL)
     (flux-deriv-replace 'c 'Ez 'EzL)
     'c
     (flux-deriv-replace 'c 'Ez 'EzL)
     'c
    '(* c c)
    (flux-deriv-replace 'By 'Ez 'EzL)
    'By
   '(* (* c c) By)
  '(* -1.0 (* (* c c) By))
  (flux-deriv-replace '(* -1.0 (* (* c c) By)) 'By 'ByL)
   (flux-deriv-replace -1.0 'By 'ByL)
   -1.0
   (flux-deriv-replace '(* (* c c) By) 'By 'ByL)
    (flux-deriv-replace '(* c c) 'By 'ByL)
     (flux-deriv-replace 'c 'By 'ByL)
     'c
     (flux-deriv-replace 'c 'By 'ByL)
     'c
    '(* c c)
    (flux-deriv-replace 'By 'By 'ByL)
    'ByL
   '(* (* c c) ByL)
  '(* -1.0 (* (* c c) ByL))
  (flux-deriv-replace '(* -1.0 (* (* c c) By)) 'Ez 'EzR)
   (flux-deriv-replace -1.0 'Ez 'EzR)
   -1.0
   (flux-deriv-replace '(* (* c c) By) 'Ez 'EzR)
    (flux-deriv-replace '(* c c) 'Ez 'EzR)
     (flux-deriv-replace 'c 'Ez 'EzR)
     'c
     (flux-deriv-replace 'c 'Ez 'EzR)
     'c
    '(* c c)
    (flux-deriv-replace 'By 'Ez 'EzR)
    'By
   '(* (* c c) By)
  '(* -1.0 (* (* c c) By))
  (flux-deriv-replace '(* -1.0 (* (* c c) By)) 'By 'ByR)
   (flux-deriv-replace -1.0 'By 'ByR)
   -1.0
   (flux-deriv-replace '(* (* c c) By) 'By 'ByR)
    (flux-deriv-replace '(* c c) 'By 'ByR)
     (flux-deriv-replace 'c 'By 'ByR)
     'c
     (flux-deriv-replace 'c 'By 'ByR)
     'c
    '(* c c)
    (flux-deriv-replace 'By 'By 'ByR)
    'ByR
   '(* (* c c) ByR)
  '(* -1.0 (* (* c c) ByR))
  (symbolic-simp '(- (* -1.0 (* (* c c) ByL)) (* -1.0 (* (* c c) ByR))))
   (symbolic-simp-rule '(- (* -1.0 (* (* c c) ByL)) (* -1.0 (* (* c c) ByR))))
   '(* -1.0 (- (* (* c c) ByL) (* (* c c) ByR)))
  (symbolic-simp '(* -1.0 (- (* (* c c) ByL) (* (* c c) ByR))))
   (symbolic-simp-rule '(* -1.0 (- (* (* c c) ByL) (* (* c c) ByR))))
    (symbolic-simp-rule -1.0)
    -1.0
    (symbolic-simp-rule '(- (* (* c c) ByL) (* (* c c) ByR)))
    '(* (* c c) (- ByL ByR))
   '(* -1.0 (* (* c c) (- ByL ByR)))
  (symbolic-simp '(* -1.0 (* (* c c) (- ByL ByR))))
   (symbolic-simp-rule '(* -1.0 (* (* c c) (- ByL ByR))))
    (symbolic-simp-rule -1.0)
    -1.0
    (symbolic-simp-rule '(* (* c c) (- ByL ByR)))
    '(* c (* c (- ByL ByR)))
   '(* -1.0 (* c (* c (- ByL ByR))))
  (symbolic-simp '(* -1.0 (* c (* c (- ByL ByR)))))
   (symbolic-simp-rule '(* -1.0 (* c (* c (- ByL ByR)))))
    (symbolic-simp-rule -1.0)
    -1.0
    (symbolic-simp-rule '(* c (* c (- ByL ByR))))
     (symbolic-simp-rule 'c)
     'c
     (symbolic-simp-rule '(* c (- ByL ByR)))
      (symbolic-simp-rule 'c)
      'c
      (symbolic-simp-rule '(- ByL ByR))
       (symbolic-simp-rule 'ByL)
       'ByL
       (symbolic-simp-rule 'ByR)
       'ByR
      '(- ByL ByR)
     '(* c (- ByL ByR))
    '(* c (* c (- ByL ByR)))
   '(* -1.0 (* c (* c (- ByL ByR))))
  '(* -1.0 (* c (* c (- ByL ByR))))
  (flux-deriv-replace '(* -1.0 Ez) 'Ez 'EzL)
   (flux-deriv-replace -1.0 'Ez 'EzL)
   -1.0
   (flux-deriv-replace 'Ez 'Ez 'EzL)
   'EzL
  '(* -1.0 EzL)
  (flux-deriv-replace '(* -1.0 EzL) 'By 'ByL)
   (flux-deriv-replace -1.0 'By 'ByL)
   -1.0
   (flux-deriv-replace 'EzL 'By 'ByL)
   'EzL
  '(* -1.0 EzL)
  (flux-deriv-replace '(* -1.0 Ez) 'Ez 'EzR)
   (flux-deriv-replace -1.0 'Ez 'EzR)
   -1.0
   (flux-deriv-replace 'Ez 'Ez 'EzR)
   'EzR
  '(* -1.0 EzR)
  (flux-deriv-replace '(* -1.0 EzR) 'By 'ByR)
   (flux-deriv-replace -1.0 'By 'ByR)
   -1.0
   (flux-deriv-replace 'EzR 'By 'ByR)
   'EzR
  '(* -1.0 EzR)
  (symbolic-simp '(- (* -1.0 EzL) (* -1.0 EzR)))
   (symbolic-simp-rule '(- (* -1.0 EzL) (* -1.0 EzR)))
   '(* -1.0 (- EzL EzR))
  (symbolic-simp '(* -1.0 (- EzL EzR)))
   (symbolic-simp-rule '(* -1.0 (- EzL EzR)))
    (symbolic-simp-rule -1.0)
    -1.0
    (symbolic-simp-rule '(- EzL EzR))
     (symbolic-simp-rule 'EzL)
     'EzL
     (symbolic-simp-rule 'EzR)
     'EzR
    '(- EzL EzR)
   '(* -1.0 (- EzL EzR))
  '(* -1.0 (- EzL EzR))
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
 #t
