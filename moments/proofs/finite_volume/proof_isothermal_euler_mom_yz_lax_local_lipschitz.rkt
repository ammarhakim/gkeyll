#lang racket

(require "../prover_core.rkt")
(require "../prover_vector.rkt")

 (prove-lax-friedrichs-vector2-1d-local-lipschitz '#hash((cons-exprs . (mom_y mom_z)) (flux-exprs . ((* mom_y u) (* mom_z u))) (max-speed-exprs . ((abs u) (abs u))) (name . "isothermal-euler-mom-yz") (parameters . ((define u 0.0)))) #:cfl 0.95 #:init-funcs '(0.0 0.0) #:nx 200 #:t-final 0.1 #:x0 0.0 #:x1 1.0)
  (symbolic-hessian '(* mom_y u) '(mom_y mom_z))
   (symbolic-gradient '(* mom_y u) '(mom_y mom_z))
    (symbolic-diff '(* mom_y u) 'mom_y)
     (symbolic-diff 'mom_y 'mom_y)
     1.0
     (symbolic-diff 'u 'mom_y)
     0.0
    '(+ (* 1.0 u) (* mom_y 0.0))
    (symbolic-simp '(+ (* 1.0 u) (* mom_y 0.0)))
     (symbolic-simp-rule '(+ (* 1.0 u) (* mom_y 0.0)))
      (symbolic-simp-rule '(* 1.0 u))
      'u
      (symbolic-simp-rule '(* mom_y 0.0))
      '(* 0.0 mom_y)
     '(+ u (* 0.0 mom_y))
    (symbolic-simp '(+ u (* 0.0 mom_y)))
     (symbolic-simp-rule '(+ u (* 0.0 mom_y)))
      (symbolic-simp-rule 'u)
      'u
      (symbolic-simp-rule '(* 0.0 mom_y))
      0.0
     '(+ u 0.0)
    (symbolic-simp '(+ u 0.0))
     (symbolic-simp-rule '(+ u 0.0))
     '(+ 0.0 u)
    (symbolic-simp '(+ 0.0 u))
     (symbolic-simp-rule '(+ 0.0 u))
     'u
    (symbolic-simp 'u)
     (symbolic-simp-rule 'u)
     'u
    'u
    (symbolic-diff '(* mom_y u) 'mom_z)
     (symbolic-diff 'mom_y 'mom_z)
     0.0
     (symbolic-diff 'u 'mom_z)
     0.0
    '(+ (* 0.0 u) (* mom_y 0.0))
    (symbolic-simp '(+ (* 0.0 u) (* mom_y 0.0)))
     (symbolic-simp-rule '(+ (* 0.0 u) (* mom_y 0.0)))
      (symbolic-simp-rule '(* 0.0 u))
      0.0
      (symbolic-simp-rule '(* mom_y 0.0))
      '(* 0.0 mom_y)
     '(+ 0.0 (* 0.0 mom_y))
    (symbolic-simp '(+ 0.0 (* 0.0 mom_y)))
     (symbolic-simp-rule '(+ 0.0 (* 0.0 mom_y)))
     '(* 0.0 mom_y)
    (symbolic-simp '(* 0.0 mom_y))
     (symbolic-simp-rule '(* 0.0 mom_y))
     0.0
    (symbolic-simp 0.0)
     (symbolic-simp-rule 0.0)
     0.0
    0.0
   '(u 0.0)
  (symbolic-jacobian '(u 0.0) '(mom_y mom_z))
   (symbolic-diff 'u 'mom_y)
   0.0
   (symbolic-simp 0.0)
    (symbolic-simp-rule 0.0)
    0.0
   0.0
   (symbolic-diff 'u 'mom_z)
   0.0
   (symbolic-simp 0.0)
    (symbolic-simp-rule 0.0)
    0.0
   0.0
   (symbolic-diff 0.0 'mom_y)
   0.0
   (symbolic-simp 0.0)
    (symbolic-simp-rule 0.0)
    0.0
   0.0
   (symbolic-diff 0.0 'mom_z)
   0.0
   (symbolic-simp 0.0)
    (symbolic-simp-rule 0.0)
    0.0
   0.0
  '((0.0 0.0) (0.0 0.0))
  (symbolic-hessian '(* mom_z u) '(mom_y mom_z))
   (symbolic-gradient '(* mom_z u) '(mom_y mom_z))
    (symbolic-diff '(* mom_z u) 'mom_y)
     (symbolic-diff 'mom_z 'mom_y)
     0.0
     (symbolic-diff 'u 'mom_y)
     0.0
    '(+ (* 0.0 u) (* mom_z 0.0))
    (symbolic-simp '(+ (* 0.0 u) (* mom_z 0.0)))
     (symbolic-simp-rule '(+ (* 0.0 u) (* mom_z 0.0)))
      (symbolic-simp-rule '(* 0.0 u))
      0.0
      (symbolic-simp-rule '(* mom_z 0.0))
      '(* 0.0 mom_z)
     '(+ 0.0 (* 0.0 mom_z))
    (symbolic-simp '(+ 0.0 (* 0.0 mom_z)))
     (symbolic-simp-rule '(+ 0.0 (* 0.0 mom_z)))
     '(* 0.0 mom_z)
    (symbolic-simp '(* 0.0 mom_z))
     (symbolic-simp-rule '(* 0.0 mom_z))
     0.0
    (symbolic-simp 0.0)
     (symbolic-simp-rule 0.0)
     0.0
    0.0
    (symbolic-diff '(* mom_z u) 'mom_z)
     (symbolic-diff 'mom_z 'mom_z)
     1.0
     (symbolic-diff 'u 'mom_z)
     0.0
    '(+ (* 1.0 u) (* mom_z 0.0))
    (symbolic-simp '(+ (* 1.0 u) (* mom_z 0.0)))
     (symbolic-simp-rule '(+ (* 1.0 u) (* mom_z 0.0)))
      (symbolic-simp-rule '(* 1.0 u))
      'u
      (symbolic-simp-rule '(* mom_z 0.0))
      '(* 0.0 mom_z)
     '(+ u (* 0.0 mom_z))
    (symbolic-simp '(+ u (* 0.0 mom_z)))
     (symbolic-simp-rule '(+ u (* 0.0 mom_z)))
      (symbolic-simp-rule 'u)
      'u
      (symbolic-simp-rule '(* 0.0 mom_z))
      0.0
     '(+ u 0.0)
    (symbolic-simp '(+ u 0.0))
     (symbolic-simp-rule '(+ u 0.0))
     '(+ 0.0 u)
    (symbolic-simp '(+ 0.0 u))
     (symbolic-simp-rule '(+ 0.0 u))
     'u
    (symbolic-simp 'u)
     (symbolic-simp-rule 'u)
     'u
    'u
   '(0.0 u)
  (symbolic-jacobian '(0.0 u) '(mom_y mom_z))
   (symbolic-diff 0.0 'mom_y)
   0.0
   (symbolic-simp 0.0)
    (symbolic-simp-rule 0.0)
    0.0
   0.0
   (symbolic-diff 0.0 'mom_z)
   0.0
   (symbolic-simp 0.0)
    (symbolic-simp-rule 0.0)
    0.0
   0.0
   (symbolic-diff 'u 'mom_y)
   0.0
   (symbolic-simp 0.0)
    (symbolic-simp-rule 0.0)
    0.0
   0.0
   (symbolic-diff 'u 'mom_z)
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
  (is-real 0.0 '(mom_y mom_z) '((define u 0.0)))
  #t
  (is-real 0.0 '(mom_y mom_z) '((define u 0.0)))
  #t
  (is-real 0.0 '(mom_y mom_z) '((define u 0.0)))
  #t
  (is-non-negative 0.0 '((define u 0.0)))
  #t
  (is-non-negative 0.0 '((define u 0.0)))
  #t
  (is-non-negative 0.0 '((define u 0.0)))
  #t
  (is-non-negative 0.0 '((define u 0.0)))
  #t
 #t
