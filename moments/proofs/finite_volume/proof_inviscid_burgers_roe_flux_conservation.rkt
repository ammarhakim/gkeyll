#lang racket

(require "../prover_core.rkt")

 (prove-roe-scalar-1d-flux-conservation '#hash((cons-expr . u) (flux-expr . (* 0.5 u u)) (max-speed-expr . (abs u)) (name . "burgers") (parameters . ())) #:cfl 0.95 #:init-func '(cond ((< (abs x) 1.0) 3.0) (else -1.0)) #:nx 200 #:t-final 0.5 #:x0 -3.0 #:x1 3.0)
  (symbolic-diff '(* 0.5 u u) 'u)
   (symbolic-diff 0.5 'u)
   0.0
   (symbolic-diff 'u 'u)
   1.0
   (symbolic-diff 'u 'u)
   1.0
  '(+ (* 0.0 u u) (* 0.5 1.0 u) (* 0.5 u 1.0))
  (symbolic-simp '(+ (* 0.0 u u) (* 0.5 1.0 u) (* 0.5 u 1.0)))
   (symbolic-simp-rule '(+ (* 0.0 u u) (* 0.5 1.0 u) (* 0.5 u 1.0)))
   '(+ (+ (* 0.0 u u) (* 0.5 1.0 u)) (* 0.5 u 1.0))
  (symbolic-simp '(+ (+ (* 0.0 u u) (* 0.5 1.0 u)) (* 0.5 u 1.0)))
   (symbolic-simp-rule '(+ (+ (* 0.0 u u) (* 0.5 1.0 u)) (* 0.5 u 1.0)))
   '(+ (* 0.0 u u) (+ (* 0.5 1.0 u) (* 0.5 u 1.0)))
  (symbolic-simp '(+ (* 0.0 u u) (+ (* 0.5 1.0 u) (* 0.5 u 1.0))))
   (symbolic-simp-rule '(+ (* 0.0 u u) (+ (* 0.5 1.0 u) (* 0.5 u 1.0))))
    (symbolic-simp-rule '(* 0.0 u u))
    '(* (* 0.0 u) u)
    (symbolic-simp-rule '(+ (* 0.5 1.0 u) (* 0.5 u 1.0)))
     (symbolic-simp-rule '(* 0.5 1.0 u))
     '(* (* 0.5 1.0) u)
     (symbolic-simp-rule '(* 0.5 u 1.0))
     '(* (* 0.5 u) 1.0)
    '(+ (* (* 0.5 1.0) u) (* (* 0.5 u) 1.0))
   '(+ (* (* 0.0 u) u) (+ (* (* 0.5 1.0) u) (* (* 0.5 u) 1.0)))
  (symbolic-simp '(+ (* (* 0.0 u) u) (+ (* (* 0.5 1.0) u) (* (* 0.5 u) 1.0))))
   (symbolic-simp-rule '(+ (* (* 0.0 u) u) (+ (* (* 0.5 1.0) u) (* (* 0.5 u) 1.0))))
    (symbolic-simp-rule '(* (* 0.0 u) u))
    '(* 0.0 (* u u))
    (symbolic-simp-rule '(+ (* (* 0.5 1.0) u) (* (* 0.5 u) 1.0)))
     (symbolic-simp-rule '(* (* 0.5 1.0) u))
     '(* 0.5 (* 1.0 u))
     (symbolic-simp-rule '(* (* 0.5 u) 1.0))
     '(* 0.5 (* u 1.0))
    '(+ (* 0.5 (* 1.0 u)) (* 0.5 (* u 1.0)))
   '(+ (* 0.0 (* u u)) (+ (* 0.5 (* 1.0 u)) (* 0.5 (* u 1.0))))
  (symbolic-simp '(+ (* 0.0 (* u u)) (+ (* 0.5 (* 1.0 u)) (* 0.5 (* u 1.0)))))
   (symbolic-simp-rule '(+ (* 0.0 (* u u)) (+ (* 0.5 (* 1.0 u)) (* 0.5 (* u 1.0)))))
    (symbolic-simp-rule '(* 0.0 (* u u)))
    0.0
    (symbolic-simp-rule '(+ (* 0.5 (* 1.0 u)) (* 0.5 (* u 1.0))))
     (symbolic-simp-rule '(* 0.5 (* 1.0 u)))
     '(* 0.5 u)
     (symbolic-simp-rule '(* 0.5 (* u 1.0)))
      (symbolic-simp-rule 0.5)
      0.5
      (symbolic-simp-rule '(* u 1.0))
      '(* 1.0 u)
     '(* 0.5 (* 1.0 u))
    '(+ (* 0.5 u) (* 0.5 (* 1.0 u)))
   '(+ 0.0 (+ (* 0.5 u) (* 0.5 (* 1.0 u))))
  (symbolic-simp '(+ 0.0 (+ (* 0.5 u) (* 0.5 (* 1.0 u)))))
   (symbolic-simp-rule '(+ 0.0 (+ (* 0.5 u) (* 0.5 (* 1.0 u)))))
   '(+ (* 0.5 u) (* 0.5 (* 1.0 u)))
  (symbolic-simp '(+ (* 0.5 u) (* 0.5 (* 1.0 u))))
   (symbolic-simp-rule '(+ (* 0.5 u) (* 0.5 (* 1.0 u))))
    (symbolic-simp-rule '(* 0.5 u))
     (symbolic-simp-rule 0.5)
     0.5
     (symbolic-simp-rule 'u)
     'u
    '(* 0.5 u)
    (symbolic-simp-rule '(* 0.5 (* 1.0 u)))
    '(* 0.5 u)
   '(+ (* 0.5 u) (* 0.5 u))
  (symbolic-simp '(+ (* 0.5 u) (* 0.5 u)))
   (symbolic-simp-rule '(+ (* 0.5 u) (* 0.5 u)))
   '(* (+ 0.5 0.5) u)
  (symbolic-simp '(* (+ 0.5 0.5) u))
   (symbolic-simp-rule '(* (+ 0.5 0.5) u))
    (symbolic-simp-rule '(+ 0.5 0.5))
    1.0
    (symbolic-simp-rule 'u)
    'u
   '(* 1.0 u)
  (symbolic-simp '(* 1.0 u))
   (symbolic-simp-rule '(* 1.0 u))
   'u
  (symbolic-simp 'u)
   (symbolic-simp-rule 'u)
   'u
  'u
  (symbolic-roe-function 'u 'u)
   (flux-deriv-replace 'u 'u 'uL)
   'uL
   (flux-deriv-replace 'u 'u 'uR)
   'uR
  (symbolic-simp '(+ (* 0.5 uL) (* 0.5 uR)))
   (symbolic-simp-rule '(+ (* 0.5 uL) (* 0.5 uR)))
    (symbolic-simp-rule '(* 0.5 uL))
     (symbolic-simp-rule 0.5)
     0.5
     (symbolic-simp-rule 'uL)
     'uL
    '(* 0.5 uL)
    (symbolic-simp-rule '(* 0.5 uR))
     (symbolic-simp-rule 0.5)
     0.5
     (symbolic-simp-rule 'uR)
     'uR
    '(* 0.5 uR)
   '(+ (* 0.5 uL) (* 0.5 uR))
  '(+ (* 0.5 uL) (* 0.5 uR))
  (symbolic-simp '(* (+ (* 0.5 uL) (* 0.5 uR)) (- uL uR)))
   (symbolic-simp-rule '(* (+ (* 0.5 uL) (* 0.5 uR)) (- uL uR)))
   '(- (* 0.5 (* uL uL)) (* 0.5 (* uR uR)))
  (symbolic-simp '(- (* 0.5 (* uL uL)) (* 0.5 (* uR uR))))
   (symbolic-simp-rule '(- (* 0.5 (* uL uL)) (* 0.5 (* uR uR))))
   '(* 0.5 (- (* uL uL) (* uR uR)))
  (symbolic-simp '(* 0.5 (- (* uL uL) (* uR uR))))
   (symbolic-simp-rule '(* 0.5 (- (* uL uL) (* uR uR))))
    (symbolic-simp-rule 0.5)
    0.5
    (symbolic-simp-rule '(- (* uL uL) (* uR uR)))
     (symbolic-simp-rule '(* uL uL))
      (symbolic-simp-rule 'uL)
      'uL
      (symbolic-simp-rule 'uL)
      'uL
     '(* uL uL)
     (symbolic-simp-rule '(* uR uR))
      (symbolic-simp-rule 'uR)
      'uR
      (symbolic-simp-rule 'uR)
      'uR
     '(* uR uR)
    '(- (* uL uL) (* uR uR))
   '(* 0.5 (- (* uL uL) (* uR uR)))
  '(* 0.5 (- (* uL uL) (* uR uR)))
  (flux-deriv-replace '(* 0.5 u u) 'u 'uL)
   (flux-deriv-replace 0.5 'u 'uL)
   0.5
   (flux-deriv-replace 'u 'u 'uL)
   'uL
   (flux-deriv-replace 'u 'u 'uL)
   'uL
  '(* 0.5 uL uL)
  (flux-deriv-replace '(* 0.5 u u) 'u 'uR)
   (flux-deriv-replace 0.5 'u 'uR)
   0.5
   (flux-deriv-replace 'u 'u 'uR)
   'uR
   (flux-deriv-replace 'u 'u 'uR)
   'uR
  '(* 0.5 uR uR)
  (symbolic-simp '(- (* 0.5 uL uL) (* 0.5 uR uR)))
   (symbolic-simp-rule '(- (* 0.5 uL uL) (* 0.5 uR uR)))
    (symbolic-simp-rule '(* 0.5 uL uL))
    '(* (* 0.5 uL) uL)
    (symbolic-simp-rule '(* 0.5 uR uR))
    '(* (* 0.5 uR) uR)
   '(- (* (* 0.5 uL) uL) (* (* 0.5 uR) uR))
  (symbolic-simp '(- (* (* 0.5 uL) uL) (* (* 0.5 uR) uR)))
   (symbolic-simp-rule '(- (* (* 0.5 uL) uL) (* (* 0.5 uR) uR)))
    (symbolic-simp-rule '(* (* 0.5 uL) uL))
    '(* 0.5 (* uL uL))
    (symbolic-simp-rule '(* (* 0.5 uR) uR))
    '(* 0.5 (* uR uR))
   '(- (* 0.5 (* uL uL)) (* 0.5 (* uR uR)))
  (symbolic-simp '(- (* 0.5 (* uL uL)) (* 0.5 (* uR uR))))
   (symbolic-simp-rule '(- (* 0.5 (* uL uL)) (* 0.5 (* uR uR))))
   '(* 0.5 (- (* uL uL) (* uR uR)))
  (symbolic-simp '(* 0.5 (- (* uL uL) (* uR uR))))
   (symbolic-simp-rule '(* 0.5 (- (* uL uL) (* uR uR))))
    (symbolic-simp-rule 0.5)
    0.5
    (symbolic-simp-rule '(- (* uL uL) (* uR uR)))
     (symbolic-simp-rule '(* uL uL))
      (symbolic-simp-rule 'uL)
      'uL
      (symbolic-simp-rule 'uL)
      'uL
     '(* uL uL)
     (symbolic-simp-rule '(* uR uR))
      (symbolic-simp-rule 'uR)
      'uR
      (symbolic-simp-rule 'uR)
      'uR
     '(* uR uR)
    '(- (* uL uL) (* uR uR))
   '(* 0.5 (- (* uL uL) (* uR uR)))
  '(* 0.5 (- (* uL uL) (* uR uR)))
  (is-real '(cond ((< (abs x) 1.0) 3.0) (else -1.0)) '(u) '())
   (is-real 3.0 '(u) '())
   #t
  (is-real -1.0 '(u) '())
  #t
 #t
