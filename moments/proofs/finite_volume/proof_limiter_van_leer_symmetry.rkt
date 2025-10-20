#lang racket

(require "../prover_core.rkt")

 (prove-flux-limiter-symmetry '#hash((limiter-expr . (/ (+ r (abs r)) (+ 1.0 (abs r)))) (limiter-ratio . r) (name . "van_leer")))
  (symbolic-simp-positive '(/ (/ (+ r (abs r)) (+ 1.0 (abs r))) r) 'r)
   (symbolic-simp-positive-rule '(/ (/ (+ r (abs r)) (+ 1.0 (abs r))) r) 'r)
    (symbolic-simp-positive-rule '(/ (+ r (abs r)) (+ 1.0 (abs r))) 'r)
     (symbolic-simp-positive-rule '(+ r (abs r)) 'r)
      (symbolic-simp-positive-rule 'r 'r)
      'r
      (symbolic-simp-positive-rule '(abs r) 'r)
      'r
     '(+ r r)
     (symbolic-simp-positive-rule '(+ 1.0 (abs r)) 'r)
      (symbolic-simp-positive-rule 1.0 'r)
      1.0
      (symbolic-simp-positive-rule '(abs r) 'r)
      'r
     '(+ 1.0 r)
    '(/ (+ r r) (+ 1.0 r))
    (symbolic-simp-positive-rule 'r 'r)
    'r
   '(/ (/ (+ r r) (+ 1.0 r)) r)
  (symbolic-simp-positive '(/ (/ (+ r r) (+ 1.0 r)) r) 'r)
   (symbolic-simp-positive-rule '(/ (/ (+ r r) (+ 1.0 r)) r) 'r)
    (symbolic-simp-positive-rule '(/ (+ r r) (+ 1.0 r)) 'r)
     (symbolic-simp-positive-rule '(+ r r) 'r)
      (symbolic-simp-positive-rule 'r 'r)
      'r
      (symbolic-simp-positive-rule 'r 'r)
      'r
     '(+ r r)
     (symbolic-simp-positive-rule '(+ 1.0 r) 'r)
      (symbolic-simp-positive-rule 1.0 'r)
      1.0
      (symbolic-simp-positive-rule 'r 'r)
      'r
     '(+ 1.0 r)
    '(/ (+ r r) (+ 1.0 r))
    (symbolic-simp-positive-rule 'r 'r)
    'r
   '(/ (/ (+ r r) (+ 1.0 r)) r)
  '(/ (/ (+ r r) (+ 1.0 r)) r)
  (symbolic-simp '(/ (/ (+ r r) (+ 1.0 r)) r))
   (symbolic-simp-rule '(/ (/ (+ r r) (+ 1.0 r)) r))
    (symbolic-simp-rule '(/ (+ r r) (+ 1.0 r)))
    '(+ (/ r (+ 1.0 r)) (/ r (+ 1.0 r)))
    (symbolic-simp-rule 'r)
    'r
   '(/ (+ (/ r (+ 1.0 r)) (/ r (+ 1.0 r))) r)
  (symbolic-simp '(/ (+ (/ r (+ 1.0 r)) (/ r (+ 1.0 r))) r))
   (symbolic-simp-rule '(/ (+ (/ r (+ 1.0 r)) (/ r (+ 1.0 r))) r))
   '(+ (/ (/ r (+ 1.0 r)) r) (/ (/ r (+ 1.0 r)) r))
  (symbolic-simp '(+ (/ (/ r (+ 1.0 r)) r) (/ (/ r (+ 1.0 r)) r)))
   (symbolic-simp-rule '(+ (/ (/ r (+ 1.0 r)) r) (/ (/ r (+ 1.0 r)) r)))
   '(* 2.0 (/ (/ r (+ 1.0 r)) r))
  (symbolic-simp '(* 2.0 (/ (/ r (+ 1.0 r)) r)))
   (symbolic-simp-rule '(* 2.0 (/ (/ r (+ 1.0 r)) r)))
    (symbolic-simp-rule 2.0)
    2.0
    (symbolic-simp-rule '(/ (/ r (+ 1.0 r)) r))
    '(/ 1.0 (+ 1.0 r))
   '(* 2.0 (/ 1.0 (+ 1.0 r)))
  (symbolic-simp '(* 2.0 (/ 1.0 (+ 1.0 r))))
   (symbolic-simp-rule '(* 2.0 (/ 1.0 (+ 1.0 r))))
   '(/ 2.0 (+ 1.0 r))
  (symbolic-simp '(/ 2.0 (+ 1.0 r)))
   (symbolic-simp-rule '(/ 2.0 (+ 1.0 r)))
    (symbolic-simp-rule 2.0)
    2.0
    (symbolic-simp-rule '(+ 1.0 r))
     (symbolic-simp-rule 1.0)
     1.0
     (symbolic-simp-rule 'r)
     'r
    '(+ 1.0 r)
   '(/ 2.0 (+ 1.0 r))
  '(/ 2.0 (+ 1.0 r))
  (symbolic-simp-positive '(/ 2.0 (+ 1.0 r)) 'r)
   (symbolic-simp-positive-rule '(/ 2.0 (+ 1.0 r)) 'r)
    (symbolic-simp-positive-rule 2.0 'r)
    2.0
    (symbolic-simp-positive-rule '(+ 1.0 r) 'r)
     (symbolic-simp-positive-rule 1.0 'r)
     1.0
     (symbolic-simp-positive-rule 'r 'r)
     'r
    '(+ 1.0 r)
   '(/ 2.0 (+ 1.0 r))
  '(/ 2.0 (+ 1.0 r))
  (symbolic-simp '(/ 2.0 (+ 1.0 r)))
   (symbolic-simp-rule '(/ 2.0 (+ 1.0 r)))
    (symbolic-simp-rule 2.0)
    2.0
    (symbolic-simp-rule '(+ 1.0 r))
     (symbolic-simp-rule 1.0)
     1.0
     (symbolic-simp-rule 'r)
     'r
    '(+ 1.0 r)
   '(/ 2.0 (+ 1.0 r))
  '(/ 2.0 (+ 1.0 r))
  (variable-transform '(/ (+ r (abs r)) (+ 1.0 (abs r))) 'r '(/ 1.0 r))
   (variable-transform '/ 'r '(/ 1.0 r))
   '/
   (variable-transform '(+ r (abs r)) 'r '(/ 1.0 r))
    (variable-transform '+ 'r '(/ 1.0 r))
    '+
    (variable-transform 'r 'r '(/ 1.0 r))
    '(/ 1.0 r)
    (variable-transform '(abs r) 'r '(/ 1.0 r))
     (variable-transform 'abs 'r '(/ 1.0 r))
     'abs
     (variable-transform 'r 'r '(/ 1.0 r))
     '(/ 1.0 r)
    '(abs (/ 1.0 r))
   '(+ (/ 1.0 r) (abs (/ 1.0 r)))
   (variable-transform '(+ 1.0 (abs r)) 'r '(/ 1.0 r))
    (variable-transform '+ 'r '(/ 1.0 r))
    '+
    (variable-transform 1.0 'r '(/ 1.0 r))
    1.0
    (variable-transform '(abs r) 'r '(/ 1.0 r))
     (variable-transform 'abs 'r '(/ 1.0 r))
     'abs
     (variable-transform 'r 'r '(/ 1.0 r))
     '(/ 1.0 r)
    '(abs (/ 1.0 r))
   '(+ 1.0 (abs (/ 1.0 r)))
  '(/ (+ (/ 1.0 r) (abs (/ 1.0 r))) (+ 1.0 (abs (/ 1.0 r))))
  (symbolic-simp-positive '(/ (+ (/ 1.0 r) (abs (/ 1.0 r))) (+ 1.0 (abs (/ 1.0 r)))) 'r)
   (symbolic-simp-positive-rule '(/ (+ (/ 1.0 r) (abs (/ 1.0 r))) (+ 1.0 (abs (/ 1.0 r)))) 'r)
    (symbolic-simp-positive-rule '(+ (/ 1.0 r) (abs (/ 1.0 r))) 'r)
     (symbolic-simp-positive-rule '(/ 1.0 r) 'r)
      (symbolic-simp-positive-rule 1.0 'r)
      1.0
      (symbolic-simp-positive-rule 'r 'r)
      'r
     '(/ 1.0 r)
     (symbolic-simp-positive-rule '(abs (/ 1.0 r)) 'r)
     '(/ 1.0 r)
    '(+ (/ 1.0 r) (/ 1.0 r))
    (symbolic-simp-positive-rule '(+ 1.0 (abs (/ 1.0 r))) 'r)
     (symbolic-simp-positive-rule 1.0 'r)
     1.0
     (symbolic-simp-positive-rule '(abs (/ 1.0 r)) 'r)
     '(/ 1.0 r)
    '(+ 1.0 (/ 1.0 r))
   '(/ (+ (/ 1.0 r) (/ 1.0 r)) (+ 1.0 (/ 1.0 r)))
  (symbolic-simp-positive '(/ (+ (/ 1.0 r) (/ 1.0 r)) (+ 1.0 (/ 1.0 r))) 'r)
   (symbolic-simp-positive-rule '(/ (+ (/ 1.0 r) (/ 1.0 r)) (+ 1.0 (/ 1.0 r))) 'r)
    (symbolic-simp-positive-rule '(+ (/ 1.0 r) (/ 1.0 r)) 'r)
     (symbolic-simp-positive-rule '(/ 1.0 r) 'r)
      (symbolic-simp-positive-rule 1.0 'r)
      1.0
      (symbolic-simp-positive-rule 'r 'r)
      'r
     '(/ 1.0 r)
     (symbolic-simp-positive-rule '(/ 1.0 r) 'r)
      (symbolic-simp-positive-rule 1.0 'r)
      1.0
      (symbolic-simp-positive-rule 'r 'r)
      'r
     '(/ 1.0 r)
    '(+ (/ 1.0 r) (/ 1.0 r))
    (symbolic-simp-positive-rule '(+ 1.0 (/ 1.0 r)) 'r)
     (symbolic-simp-positive-rule 1.0 'r)
     1.0
     (symbolic-simp-positive-rule '(/ 1.0 r) 'r)
      (symbolic-simp-positive-rule 1.0 'r)
      1.0
      (symbolic-simp-positive-rule 'r 'r)
      'r
     '(/ 1.0 r)
    '(+ 1.0 (/ 1.0 r))
   '(/ (+ (/ 1.0 r) (/ 1.0 r)) (+ 1.0 (/ 1.0 r)))
  '(/ (+ (/ 1.0 r) (/ 1.0 r)) (+ 1.0 (/ 1.0 r)))
  (symbolic-simp '(/ (+ (/ 1.0 r) (/ 1.0 r)) (+ 1.0 (/ 1.0 r))))
   (symbolic-simp-rule '(/ (+ (/ 1.0 r) (/ 1.0 r)) (+ 1.0 (/ 1.0 r))))
   '(+ (/ (/ 1.0 r) (+ 1.0 (/ 1.0 r))) (/ (/ 1.0 r) (+ 1.0 (/ 1.0 r))))
  (symbolic-simp '(+ (/ (/ 1.0 r) (+ 1.0 (/ 1.0 r))) (/ (/ 1.0 r) (+ 1.0 (/ 1.0 r)))))
   (symbolic-simp-rule '(+ (/ (/ 1.0 r) (+ 1.0 (/ 1.0 r))) (/ (/ 1.0 r) (+ 1.0 (/ 1.0 r)))))
   '(* 2.0 (/ (/ 1.0 r) (+ 1.0 (/ 1.0 r))))
  (symbolic-simp '(* 2.0 (/ (/ 1.0 r) (+ 1.0 (/ 1.0 r)))))
   (symbolic-simp-rule '(* 2.0 (/ (/ 1.0 r) (+ 1.0 (/ 1.0 r)))))
    (symbolic-simp-rule 2.0)
    2.0
    (symbolic-simp-rule '(/ (/ 1.0 r) (+ 1.0 (/ 1.0 r))))
    '(/ 1.0 (+ (* 1.0 r) 1.0))
   '(* 2.0 (/ 1.0 (+ (* 1.0 r) 1.0)))
  (symbolic-simp '(* 2.0 (/ 1.0 (+ (* 1.0 r) 1.0))))
   (symbolic-simp-rule '(* 2.0 (/ 1.0 (+ (* 1.0 r) 1.0))))
   '(/ 2.0 (+ (* 1.0 r) 1.0))
  (symbolic-simp '(/ 2.0 (+ (* 1.0 r) 1.0)))
   (symbolic-simp-rule '(/ 2.0 (+ (* 1.0 r) 1.0)))
    (symbolic-simp-rule 2.0)
    2.0
    (symbolic-simp-rule '(+ (* 1.0 r) 1.0))
    '(+ 1.0 (* 1.0 r))
   '(/ 2.0 (+ 1.0 (* 1.0 r)))
  (symbolic-simp '(/ 2.0 (+ 1.0 (* 1.0 r))))
   (symbolic-simp-rule '(/ 2.0 (+ 1.0 (* 1.0 r))))
    (symbolic-simp-rule 2.0)
    2.0
    (symbolic-simp-rule '(+ 1.0 (* 1.0 r)))
     (symbolic-simp-rule 1.0)
     1.0
     (symbolic-simp-rule '(* 1.0 r))
     'r
    '(+ 1.0 r)
   '(/ 2.0 (+ 1.0 r))
  (symbolic-simp '(/ 2.0 (+ 1.0 r)))
   (symbolic-simp-rule '(/ 2.0 (+ 1.0 r)))
    (symbolic-simp-rule 2.0)
    2.0
    (symbolic-simp-rule '(+ 1.0 r))
     (symbolic-simp-rule 1.0)
     1.0
     (symbolic-simp-rule 'r)
     'r
    '(+ 1.0 r)
   '(/ 2.0 (+ 1.0 r))
  '(/ 2.0 (+ 1.0 r))
  (symbolic-simp-positive '(/ 2.0 (+ 1.0 r)) 'r)
   (symbolic-simp-positive-rule '(/ 2.0 (+ 1.0 r)) 'r)
    (symbolic-simp-positive-rule 2.0 'r)
    2.0
    (symbolic-simp-positive-rule '(+ 1.0 r) 'r)
     (symbolic-simp-positive-rule 1.0 'r)
     1.0
     (symbolic-simp-positive-rule 'r 'r)
     'r
    '(+ 1.0 r)
   '(/ 2.0 (+ 1.0 r))
  '(/ 2.0 (+ 1.0 r))
  (symbolic-simp '(/ 2.0 (+ 1.0 r)))
   (symbolic-simp-rule '(/ 2.0 (+ 1.0 r)))
    (symbolic-simp-rule 2.0)
    2.0
    (symbolic-simp-rule '(+ 1.0 r))
     (symbolic-simp-rule 1.0)
     1.0
     (symbolic-simp-rule 'r)
     'r
    '(+ 1.0 r)
   '(/ 2.0 (+ 1.0 r))
  '(/ 2.0 (+ 1.0 r))
 #t
