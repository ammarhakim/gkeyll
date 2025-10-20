#lang racket

(require "../prover_core.rkt")

 (prove-roe-scalar-1d-flux-conservation '#hash((cons-expr . u) (flux-expr . (* a u)) (max-speed-expr . (abs a)) (name . "advect") (parameters . ((define a 1.0)))) #:cfl 0.95 #:init-func '(cond ((< x 1.0) 1.0) (else 0.0)) #:nx 200 #:t-final 0.5 #:x0 0.0 #:x1 2.0)
  (symbolic-diff '(* a u) 'u)
   (symbolic-diff 'a 'u)
   0.0
   (symbolic-diff 'u 'u)
   1.0
  '(+ (* 0.0 u) (* a 1.0))
  (symbolic-simp '(+ (* 0.0 u) (* a 1.0)))
   (symbolic-simp-rule '(+ (* 0.0 u) (* a 1.0)))
    (symbolic-simp-rule '(* 0.0 u))
    0.0
    (symbolic-simp-rule '(* a 1.0))
    '(* 1.0 a)
   '(+ 0.0 (* 1.0 a))
  (symbolic-simp '(+ 0.0 (* 1.0 a)))
   (symbolic-simp-rule '(+ 0.0 (* 1.0 a)))
   '(* 1.0 a)
  (symbolic-simp '(* 1.0 a))
   (symbolic-simp-rule '(* 1.0 a))
   'a
  (symbolic-simp 'a)
   (symbolic-simp-rule 'a)
   'a
  'a
  (symbolic-roe-function 'a 'u)
   (flux-deriv-replace 'a 'u 'uL)
   'a
   (flux-deriv-replace 'a 'u 'uR)
   'a
  (symbolic-simp '(+ (* 0.5 a) (* 0.5 a)))
   (symbolic-simp-rule '(+ (* 0.5 a) (* 0.5 a)))
   '(* (+ 0.5 0.5) a)
  (symbolic-simp '(* (+ 0.5 0.5) a))
   (symbolic-simp-rule '(* (+ 0.5 0.5) a))
    (symbolic-simp-rule '(+ 0.5 0.5))
    1.0
    (symbolic-simp-rule 'a)
    'a
   '(* 1.0 a)
  (symbolic-simp '(* 1.0 a))
   (symbolic-simp-rule '(* 1.0 a))
   'a
  (symbolic-simp 'a)
   (symbolic-simp-rule 'a)
   'a
  'a
  (symbolic-simp '(* a (- uL uR)))
   (symbolic-simp-rule '(* a (- uL uR)))
    (symbolic-simp-rule 'a)
    'a
    (symbolic-simp-rule '(- uL uR))
     (symbolic-simp-rule 'uL)
     'uL
     (symbolic-simp-rule 'uR)
     'uR
    '(- uL uR)
   '(* a (- uL uR))
  '(* a (- uL uR))
  (flux-deriv-replace '(* a u) 'u 'uL)
   (flux-deriv-replace 'a 'u 'uL)
   'a
   (flux-deriv-replace 'u 'u 'uL)
   'uL
  '(* a uL)
  (flux-deriv-replace '(* a u) 'u 'uR)
   (flux-deriv-replace 'a 'u 'uR)
   'a
   (flux-deriv-replace 'u 'u 'uR)
   'uR
  '(* a uR)
  (symbolic-simp '(- (* a uL) (* a uR)))
   (symbolic-simp-rule '(- (* a uL) (* a uR)))
   '(* a (- uL uR))
  (symbolic-simp '(* a (- uL uR)))
   (symbolic-simp-rule '(* a (- uL uR)))
    (symbolic-simp-rule 'a)
    'a
    (symbolic-simp-rule '(- uL uR))
     (symbolic-simp-rule 'uL)
     'uL
     (symbolic-simp-rule 'uR)
     'uR
    '(- uL uR)
   '(* a (- uL uR))
  '(* a (- uL uR))
  (is-real 1.0 '(u) '((define a 1.0)))
  #t
  (is-real '(cond ((< x 1.0) 1.0) (else 0.0)) '(u) '((define a 1.0)))
   (is-real 1.0 '(u) '((define a 1.0)))
   #t
  (is-real 0.0 '(u) '((define a 1.0)))
  #t
 #t
