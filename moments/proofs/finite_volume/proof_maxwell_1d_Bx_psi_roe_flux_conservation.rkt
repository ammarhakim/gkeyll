#lang racket

(require "../prover_core.rkt")
(require "../prover_vector.rkt")

 (prove-roe-vector2-1d-flux-conservation '#hash((cons-exprs . (Bx psi)) (flux-exprs . ((* b_fact psi) (* b_fact (* (* c c) Bx)))) (max-speed-exprs . ((abs (* b_fact c)) (abs (* b_fact c)))) (name . "maxwell_Bxpsi") (parameters . ((define c 1.0) (define e_fact 1.0) (define b_fact 1.0)))) #:cfl 0.95 #:init-funcs '(0.0 (cond ((< x 0.0) 0.5) (else -0.5))) #:nx 200 #:t-final 1.0 #:x0 -1.5 #:x1 1.5)
  (symbolic-jacobian '((* b_fact psi) (* b_fact (* (* c c) Bx))) '(Bx psi))
   (symbolic-diff '(* b_fact psi) 'Bx)
    (symbolic-diff 'b_fact 'Bx)
    0.0
    (symbolic-diff 'psi 'Bx)
    0.0
   '(+ (* 0.0 psi) (* b_fact 0.0))
   (symbolic-simp '(+ (* 0.0 psi) (* b_fact 0.0)))
    (symbolic-simp-rule '(+ (* 0.0 psi) (* b_fact 0.0)))
     (symbolic-simp-rule '(* 0.0 psi))
     0.0
     (symbolic-simp-rule '(* b_fact 0.0))
     '(* 0.0 b_fact)
    '(+ 0.0 (* 0.0 b_fact))
   (symbolic-simp '(+ 0.0 (* 0.0 b_fact)))
    (symbolic-simp-rule '(+ 0.0 (* 0.0 b_fact)))
    '(* 0.0 b_fact)
   (symbolic-simp '(* 0.0 b_fact))
    (symbolic-simp-rule '(* 0.0 b_fact))
    0.0
   (symbolic-simp 0.0)
    (symbolic-simp-rule 0.0)
    0.0
   0.0
   (symbolic-diff '(* b_fact psi) 'psi)
    (symbolic-diff 'b_fact 'psi)
    0.0
    (symbolic-diff 'psi 'psi)
    1.0
   '(+ (* 0.0 psi) (* b_fact 1.0))
   (symbolic-simp '(+ (* 0.0 psi) (* b_fact 1.0)))
    (symbolic-simp-rule '(+ (* 0.0 psi) (* b_fact 1.0)))
     (symbolic-simp-rule '(* 0.0 psi))
     0.0
     (symbolic-simp-rule '(* b_fact 1.0))
     '(* 1.0 b_fact)
    '(+ 0.0 (* 1.0 b_fact))
   (symbolic-simp '(+ 0.0 (* 1.0 b_fact)))
    (symbolic-simp-rule '(+ 0.0 (* 1.0 b_fact)))
    '(* 1.0 b_fact)
   (symbolic-simp '(* 1.0 b_fact))
    (symbolic-simp-rule '(* 1.0 b_fact))
    'b_fact
   (symbolic-simp 'b_fact)
    (symbolic-simp-rule 'b_fact)
    'b_fact
   'b_fact
   (symbolic-diff '(* b_fact (* (* c c) Bx)) 'Bx)
    (symbolic-diff 'b_fact 'Bx)
    0.0
    (symbolic-diff '(* (* c c) Bx) 'Bx)
     (symbolic-diff '(* c c) 'Bx)
      (symbolic-diff 'c 'Bx)
      0.0
      (symbolic-diff 'c 'Bx)
      0.0
     '(+ (* 0.0 c) (* c 0.0))
     (symbolic-diff 'Bx 'Bx)
     1.0
    '(+ (* (+ (* 0.0 c) (* c 0.0)) Bx) (* (* c c) 1.0))
   '(+ (* 0.0 (* (* c c) Bx)) (* b_fact (+ (* (+ (* 0.0 c) (* c 0.0)) Bx) (* (* c c) 1.0))))
   (symbolic-simp '(+ (* 0.0 (* (* c c) Bx)) (* b_fact (+ (* (+ (* 0.0 c) (* c 0.0)) Bx) (* (* c c) 1.0)))))
    (symbolic-simp-rule '(+ (* 0.0 (* (* c c) Bx)) (* b_fact (+ (* (+ (* 0.0 c) (* c 0.0)) Bx) (* (* c c) 1.0)))))
     (symbolic-simp-rule '(* 0.0 (* (* c c) Bx)))
     0.0
     (symbolic-simp-rule '(* b_fact (+ (* (+ (* 0.0 c) (* c 0.0)) Bx) (* (* c c) 1.0))))
     '(+ (* b_fact (* (+ (* 0.0 c) (* c 0.0)) Bx)) (* b_fact (* (* c c) 1.0)))
    '(+ 0.0 (+ (* b_fact (* (+ (* 0.0 c) (* c 0.0)) Bx)) (* b_fact (* (* c c) 1.0))))
   (symbolic-simp '(+ 0.0 (+ (* b_fact (* (+ (* 0.0 c) (* c 0.0)) Bx)) (* b_fact (* (* c c) 1.0)))))
    (symbolic-simp-rule '(+ 0.0 (+ (* b_fact (* (+ (* 0.0 c) (* c 0.0)) Bx)) (* b_fact (* (* c c) 1.0)))))
    '(+ (* b_fact (* (+ (* 0.0 c) (* c 0.0)) Bx)) (* b_fact (* (* c c) 1.0)))
   (symbolic-simp '(+ (* b_fact (* (+ (* 0.0 c) (* c 0.0)) Bx)) (* b_fact (* (* c c) 1.0))))
    (symbolic-simp-rule '(+ (* b_fact (* (+ (* 0.0 c) (* c 0.0)) Bx)) (* b_fact (* (* c c) 1.0))))
     (symbolic-simp-rule '(* b_fact (* (+ (* 0.0 c) (* c 0.0)) Bx)))
      (symbolic-simp-rule 'b_fact)
      'b_fact
      (symbolic-simp-rule '(* (+ (* 0.0 c) (* c 0.0)) Bx))
       (symbolic-simp-rule '(+ (* 0.0 c) (* c 0.0)))
        (symbolic-simp-rule '(* 0.0 c))
        0.0
        (symbolic-simp-rule '(* c 0.0))
        '(* 0.0 c)
       '(+ 0.0 (* 0.0 c))
       (symbolic-simp-rule 'Bx)
       'Bx
      '(* (+ 0.0 (* 0.0 c)) Bx)
     '(* b_fact (* (+ 0.0 (* 0.0 c)) Bx))
     (symbolic-simp-rule '(* b_fact (* (* c c) 1.0)))
      (symbolic-simp-rule 'b_fact)
      'b_fact
      (symbolic-simp-rule '(* (* c c) 1.0))
      '(* c (* c 1.0))
     '(* b_fact (* c (* c 1.0)))
    '(+ (* b_fact (* (+ 0.0 (* 0.0 c)) Bx)) (* b_fact (* c (* c 1.0))))
   (symbolic-simp '(+ (* b_fact (* (+ 0.0 (* 0.0 c)) Bx)) (* b_fact (* c (* c 1.0)))))
    (symbolic-simp-rule '(+ (* b_fact (* (+ 0.0 (* 0.0 c)) Bx)) (* b_fact (* c (* c 1.0)))))
     (symbolic-simp-rule '(* b_fact (* (+ 0.0 (* 0.0 c)) Bx)))
      (symbolic-simp-rule 'b_fact)
      'b_fact
      (symbolic-simp-rule '(* (+ 0.0 (* 0.0 c)) Bx))
       (symbolic-simp-rule '(+ 0.0 (* 0.0 c)))
       '(* 0.0 c)
       (symbolic-simp-rule 'Bx)
       'Bx
      '(* (* 0.0 c) Bx)
     '(* b_fact (* (* 0.0 c) Bx))
     (symbolic-simp-rule '(* b_fact (* c (* c 1.0))))
      (symbolic-simp-rule 'b_fact)
      'b_fact
      (symbolic-simp-rule '(* c (* c 1.0)))
       (symbolic-simp-rule 'c)
       'c
       (symbolic-simp-rule '(* c 1.0))
       '(* 1.0 c)
      '(* c (* 1.0 c))
     '(* b_fact (* c (* 1.0 c)))
    '(+ (* b_fact (* (* 0.0 c) Bx)) (* b_fact (* c (* 1.0 c))))
   (symbolic-simp '(+ (* b_fact (* (* 0.0 c) Bx)) (* b_fact (* c (* 1.0 c)))))
    (symbolic-simp-rule '(+ (* b_fact (* (* 0.0 c) Bx)) (* b_fact (* c (* 1.0 c)))))
     (symbolic-simp-rule '(* b_fact (* (* 0.0 c) Bx)))
      (symbolic-simp-rule 'b_fact)
      'b_fact
      (symbolic-simp-rule '(* (* 0.0 c) Bx))
      '(* 0.0 (* c Bx))
     '(* b_fact (* 0.0 (* c Bx)))
     (symbolic-simp-rule '(* b_fact (* c (* 1.0 c))))
      (symbolic-simp-rule 'b_fact)
      'b_fact
      (symbolic-simp-rule '(* c (* 1.0 c)))
      '(* 1.0 (* c c))
     '(* b_fact (* 1.0 (* c c)))
    '(+ (* b_fact (* 0.0 (* c Bx))) (* b_fact (* 1.0 (* c c))))
   (symbolic-simp '(+ (* b_fact (* 0.0 (* c Bx))) (* b_fact (* 1.0 (* c c)))))
    (symbolic-simp-rule '(+ (* b_fact (* 0.0 (* c Bx))) (* b_fact (* 1.0 (* c c)))))
     (symbolic-simp-rule '(* b_fact (* 0.0 (* c Bx))))
     '(* 0.0 (* b_fact (* c Bx)))
     (symbolic-simp-rule '(* b_fact (* 1.0 (* c c))))
     '(* 1.0 (* b_fact (* c c)))
    '(+ (* 0.0 (* b_fact (* c Bx))) (* 1.0 (* b_fact (* c c))))
   (symbolic-simp '(+ (* 0.0 (* b_fact (* c Bx))) (* 1.0 (* b_fact (* c c)))))
    (symbolic-simp-rule '(+ (* 0.0 (* b_fact (* c Bx))) (* 1.0 (* b_fact (* c c)))))
     (symbolic-simp-rule '(* 0.0 (* b_fact (* c Bx))))
     0.0
     (symbolic-simp-rule '(* 1.0 (* b_fact (* c c))))
     '(* b_fact (* c c))
    '(+ 0.0 (* b_fact (* c c)))
   (symbolic-simp '(+ 0.0 (* b_fact (* c c))))
    (symbolic-simp-rule '(+ 0.0 (* b_fact (* c c))))
    '(* b_fact (* c c))
   (symbolic-simp '(* b_fact (* c c)))
    (symbolic-simp-rule '(* b_fact (* c c)))
     (symbolic-simp-rule 'b_fact)
     'b_fact
     (symbolic-simp-rule '(* c c))
      (symbolic-simp-rule 'c)
      'c
      (symbolic-simp-rule 'c)
      'c
     '(* c c)
    '(* b_fact (* c c))
   '(* b_fact (* c c))
   (symbolic-diff '(* b_fact (* (* c c) Bx)) 'psi)
    (symbolic-diff 'b_fact 'psi)
    0.0
    (symbolic-diff '(* (* c c) Bx) 'psi)
     (symbolic-diff '(* c c) 'psi)
      (symbolic-diff 'c 'psi)
      0.0
      (symbolic-diff 'c 'psi)
      0.0
     '(+ (* 0.0 c) (* c 0.0))
     (symbolic-diff 'Bx 'psi)
     0.0
    '(+ (* (+ (* 0.0 c) (* c 0.0)) Bx) (* (* c c) 0.0))
   '(+ (* 0.0 (* (* c c) Bx)) (* b_fact (+ (* (+ (* 0.0 c) (* c 0.0)) Bx) (* (* c c) 0.0))))
   (symbolic-simp '(+ (* 0.0 (* (* c c) Bx)) (* b_fact (+ (* (+ (* 0.0 c) (* c 0.0)) Bx) (* (* c c) 0.0)))))
    (symbolic-simp-rule '(+ (* 0.0 (* (* c c) Bx)) (* b_fact (+ (* (+ (* 0.0 c) (* c 0.0)) Bx) (* (* c c) 0.0)))))
     (symbolic-simp-rule '(* 0.0 (* (* c c) Bx)))
     0.0
     (symbolic-simp-rule '(* b_fact (+ (* (+ (* 0.0 c) (* c 0.0)) Bx) (* (* c c) 0.0))))
     '(+ (* b_fact (* (+ (* 0.0 c) (* c 0.0)) Bx)) (* b_fact (* (* c c) 0.0)))
    '(+ 0.0 (+ (* b_fact (* (+ (* 0.0 c) (* c 0.0)) Bx)) (* b_fact (* (* c c) 0.0))))
   (symbolic-simp '(+ 0.0 (+ (* b_fact (* (+ (* 0.0 c) (* c 0.0)) Bx)) (* b_fact (* (* c c) 0.0)))))
    (symbolic-simp-rule '(+ 0.0 (+ (* b_fact (* (+ (* 0.0 c) (* c 0.0)) Bx)) (* b_fact (* (* c c) 0.0)))))
    '(+ (* b_fact (* (+ (* 0.0 c) (* c 0.0)) Bx)) (* b_fact (* (* c c) 0.0)))
   (symbolic-simp '(+ (* b_fact (* (+ (* 0.0 c) (* c 0.0)) Bx)) (* b_fact (* (* c c) 0.0))))
    (symbolic-simp-rule '(+ (* b_fact (* (+ (* 0.0 c) (* c 0.0)) Bx)) (* b_fact (* (* c c) 0.0))))
     (symbolic-simp-rule '(* b_fact (* (+ (* 0.0 c) (* c 0.0)) Bx)))
      (symbolic-simp-rule 'b_fact)
      'b_fact
      (symbolic-simp-rule '(* (+ (* 0.0 c) (* c 0.0)) Bx))
       (symbolic-simp-rule '(+ (* 0.0 c) (* c 0.0)))
        (symbolic-simp-rule '(* 0.0 c))
        0.0
        (symbolic-simp-rule '(* c 0.0))
        '(* 0.0 c)
       '(+ 0.0 (* 0.0 c))
       (symbolic-simp-rule 'Bx)
       'Bx
      '(* (+ 0.0 (* 0.0 c)) Bx)
     '(* b_fact (* (+ 0.0 (* 0.0 c)) Bx))
     (symbolic-simp-rule '(* b_fact (* (* c c) 0.0)))
      (symbolic-simp-rule 'b_fact)
      'b_fact
      (symbolic-simp-rule '(* (* c c) 0.0))
      '(* c (* c 0.0))
     '(* b_fact (* c (* c 0.0)))
    '(+ (* b_fact (* (+ 0.0 (* 0.0 c)) Bx)) (* b_fact (* c (* c 0.0))))
   (symbolic-simp '(+ (* b_fact (* (+ 0.0 (* 0.0 c)) Bx)) (* b_fact (* c (* c 0.0)))))
    (symbolic-simp-rule '(+ (* b_fact (* (+ 0.0 (* 0.0 c)) Bx)) (* b_fact (* c (* c 0.0)))))
     (symbolic-simp-rule '(* b_fact (* (+ 0.0 (* 0.0 c)) Bx)))
      (symbolic-simp-rule 'b_fact)
      'b_fact
      (symbolic-simp-rule '(* (+ 0.0 (* 0.0 c)) Bx))
       (symbolic-simp-rule '(+ 0.0 (* 0.0 c)))
       '(* 0.0 c)
       (symbolic-simp-rule 'Bx)
       'Bx
      '(* (* 0.0 c) Bx)
     '(* b_fact (* (* 0.0 c) Bx))
     (symbolic-simp-rule '(* b_fact (* c (* c 0.0))))
      (symbolic-simp-rule 'b_fact)
      'b_fact
      (symbolic-simp-rule '(* c (* c 0.0)))
       (symbolic-simp-rule 'c)
       'c
       (symbolic-simp-rule '(* c 0.0))
       '(* 0.0 c)
      '(* c (* 0.0 c))
     '(* b_fact (* c (* 0.0 c)))
    '(+ (* b_fact (* (* 0.0 c) Bx)) (* b_fact (* c (* 0.0 c))))
   (symbolic-simp '(+ (* b_fact (* (* 0.0 c) Bx)) (* b_fact (* c (* 0.0 c)))))
    (symbolic-simp-rule '(+ (* b_fact (* (* 0.0 c) Bx)) (* b_fact (* c (* 0.0 c)))))
     (symbolic-simp-rule '(* b_fact (* (* 0.0 c) Bx)))
      (symbolic-simp-rule 'b_fact)
      'b_fact
      (symbolic-simp-rule '(* (* 0.0 c) Bx))
      '(* 0.0 (* c Bx))
     '(* b_fact (* 0.0 (* c Bx)))
     (symbolic-simp-rule '(* b_fact (* c (* 0.0 c))))
      (symbolic-simp-rule 'b_fact)
      'b_fact
      (symbolic-simp-rule '(* c (* 0.0 c)))
      '(* 0.0 (* c c))
     '(* b_fact (* 0.0 (* c c)))
    '(+ (* b_fact (* 0.0 (* c Bx))) (* b_fact (* 0.0 (* c c))))
   (symbolic-simp '(+ (* b_fact (* 0.0 (* c Bx))) (* b_fact (* 0.0 (* c c)))))
    (symbolic-simp-rule '(+ (* b_fact (* 0.0 (* c Bx))) (* b_fact (* 0.0 (* c c)))))
     (symbolic-simp-rule '(* b_fact (* 0.0 (* c Bx))))
     '(* 0.0 (* b_fact (* c Bx)))
     (symbolic-simp-rule '(* b_fact (* 0.0 (* c c))))
     '(* 0.0 (* b_fact (* c c)))
    '(+ (* 0.0 (* b_fact (* c Bx))) (* 0.0 (* b_fact (* c c))))
   (symbolic-simp '(+ (* 0.0 (* b_fact (* c Bx))) (* 0.0 (* b_fact (* c c)))))
    (symbolic-simp-rule '(+ (* 0.0 (* b_fact (* c Bx))) (* 0.0 (* b_fact (* c c)))))
     (symbolic-simp-rule '(* 0.0 (* b_fact (* c Bx))))
     0.0
     (symbolic-simp-rule '(* 0.0 (* b_fact (* c c))))
     0.0
    '(+ 0.0 0.0)
   (symbolic-simp '(+ 0.0 0.0))
    (symbolic-simp-rule '(+ 0.0 0.0))
    0.0
   (symbolic-simp 0.0)
    (symbolic-simp-rule 0.0)
    0.0
   0.0
  '((0.0 b_fact) ((* b_fact (* c c)) 0.0))
  (symbolic-roe-matrix '((0.0 b_fact) ((* b_fact (* c c)) 0.0)) '(Bx psi))
   (flux-deriv-replace 0.0 'Bx 'BxL)
   0.0
   (flux-deriv-replace 0.0 'psi 'psiL)
   0.0
   (flux-deriv-replace 0.0 'Bx 'BxR)
   0.0
   (flux-deriv-replace 0.0 'psi 'psiR)
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
   (flux-deriv-replace 'b_fact 'Bx 'BxL)
   'b_fact
   (flux-deriv-replace 'b_fact 'psi 'psiL)
   'b_fact
   (flux-deriv-replace 'b_fact 'Bx 'BxR)
   'b_fact
   (flux-deriv-replace 'b_fact 'psi 'psiR)
   'b_fact
   (symbolic-simp '(+ (* 0.5 b_fact) (* 0.5 b_fact)))
    (symbolic-simp-rule '(+ (* 0.5 b_fact) (* 0.5 b_fact)))
    '(* (+ 0.5 0.5) b_fact)
   (symbolic-simp '(* (+ 0.5 0.5) b_fact))
    (symbolic-simp-rule '(* (+ 0.5 0.5) b_fact))
     (symbolic-simp-rule '(+ 0.5 0.5))
     1.0
     (symbolic-simp-rule 'b_fact)
     'b_fact
    '(* 1.0 b_fact)
   (symbolic-simp '(* 1.0 b_fact))
    (symbolic-simp-rule '(* 1.0 b_fact))
    'b_fact
   (symbolic-simp 'b_fact)
    (symbolic-simp-rule 'b_fact)
    'b_fact
   'b_fact
   (flux-deriv-replace '(* b_fact (* c c)) 'Bx 'BxL)
    (flux-deriv-replace 'b_fact 'Bx 'BxL)
    'b_fact
    (flux-deriv-replace '(* c c) 'Bx 'BxL)
     (flux-deriv-replace 'c 'Bx 'BxL)
     'c
     (flux-deriv-replace 'c 'Bx 'BxL)
     'c
    '(* c c)
   '(* b_fact (* c c))
   (flux-deriv-replace '(* b_fact (* c c)) 'psi 'psiL)
    (flux-deriv-replace 'b_fact 'psi 'psiL)
    'b_fact
    (flux-deriv-replace '(* c c) 'psi 'psiL)
     (flux-deriv-replace 'c 'psi 'psiL)
     'c
     (flux-deriv-replace 'c 'psi 'psiL)
     'c
    '(* c c)
   '(* b_fact (* c c))
   (flux-deriv-replace '(* b_fact (* c c)) 'Bx 'BxR)
    (flux-deriv-replace 'b_fact 'Bx 'BxR)
    'b_fact
    (flux-deriv-replace '(* c c) 'Bx 'BxR)
     (flux-deriv-replace 'c 'Bx 'BxR)
     'c
     (flux-deriv-replace 'c 'Bx 'BxR)
     'c
    '(* c c)
   '(* b_fact (* c c))
   (flux-deriv-replace '(* b_fact (* c c)) 'psi 'psiR)
    (flux-deriv-replace 'b_fact 'psi 'psiR)
    'b_fact
    (flux-deriv-replace '(* c c) 'psi 'psiR)
     (flux-deriv-replace 'c 'psi 'psiR)
     'c
     (flux-deriv-replace 'c 'psi 'psiR)
     'c
    '(* c c)
   '(* b_fact (* c c))
   (symbolic-simp '(+ (* 0.5 (* b_fact (* c c))) (* 0.5 (* b_fact (* c c)))))
    (symbolic-simp-rule '(+ (* 0.5 (* b_fact (* c c))) (* 0.5 (* b_fact (* c c)))))
    '(* (+ 0.5 0.5) (* b_fact (* c c)))
   (symbolic-simp '(* (+ 0.5 0.5) (* b_fact (* c c))))
    (symbolic-simp-rule '(* (+ 0.5 0.5) (* b_fact (* c c))))
     (symbolic-simp-rule '(+ 0.5 0.5))
     1.0
     (symbolic-simp-rule '(* b_fact (* c c)))
      (symbolic-simp-rule 'b_fact)
      'b_fact
      (symbolic-simp-rule '(* c c))
       (symbolic-simp-rule 'c)
       'c
       (symbolic-simp-rule 'c)
       'c
      '(* c c)
     '(* b_fact (* c c))
    '(* 1.0 (* b_fact (* c c)))
   (symbolic-simp '(* 1.0 (* b_fact (* c c))))
    (symbolic-simp-rule '(* 1.0 (* b_fact (* c c))))
    '(* b_fact (* c c))
   (symbolic-simp '(* b_fact (* c c)))
    (symbolic-simp-rule '(* b_fact (* c c)))
     (symbolic-simp-rule 'b_fact)
     'b_fact
     (symbolic-simp-rule '(* c c))
      (symbolic-simp-rule 'c)
      'c
      (symbolic-simp-rule 'c)
      'c
     '(* c c)
    '(* b_fact (* c c))
   '(* b_fact (* c c))
   (flux-deriv-replace 0.0 'Bx 'BxL)
   0.0
   (flux-deriv-replace 0.0 'psi 'psiL)
   0.0
   (flux-deriv-replace 0.0 'Bx 'BxR)
   0.0
   (flux-deriv-replace 0.0 'psi 'psiR)
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
  '((0.0 b_fact) ((* b_fact (* c c)) 0.0))
  (symbolic-simp '(- BxL BxR))
   (symbolic-simp-rule '(- BxL BxR))
    (symbolic-simp-rule 'BxL)
    'BxL
    (symbolic-simp-rule 'BxR)
    'BxR
   '(- BxL BxR)
  '(- BxL BxR)
  (symbolic-simp '(- psiL psiR))
   (symbolic-simp-rule '(- psiL psiR))
    (symbolic-simp-rule 'psiL)
    'psiL
    (symbolic-simp-rule 'psiR)
    'psiR
   '(- psiL psiR)
  '(- psiL psiR)
  (symbolic-simp '(+ (* 0.0 (- BxL BxR)) (* b_fact (- psiL psiR))))
   (symbolic-simp-rule '(+ (* 0.0 (- BxL BxR)) (* b_fact (- psiL psiR))))
    (symbolic-simp-rule '(* 0.0 (- BxL BxR)))
    0.0
    (symbolic-simp-rule '(* b_fact (- psiL psiR)))
     (symbolic-simp-rule 'b_fact)
     'b_fact
     (symbolic-simp-rule '(- psiL psiR))
      (symbolic-simp-rule 'psiL)
      'psiL
      (symbolic-simp-rule 'psiR)
      'psiR
     '(- psiL psiR)
    '(* b_fact (- psiL psiR))
   '(+ 0.0 (* b_fact (- psiL psiR)))
  (symbolic-simp '(+ 0.0 (* b_fact (- psiL psiR))))
   (symbolic-simp-rule '(+ 0.0 (* b_fact (- psiL psiR))))
   '(* b_fact (- psiL psiR))
  (symbolic-simp '(* b_fact (- psiL psiR)))
   (symbolic-simp-rule '(* b_fact (- psiL psiR)))
    (symbolic-simp-rule 'b_fact)
    'b_fact
    (symbolic-simp-rule '(- psiL psiR))
     (symbolic-simp-rule 'psiL)
     'psiL
     (symbolic-simp-rule 'psiR)
     'psiR
    '(- psiL psiR)
   '(* b_fact (- psiL psiR))
  '(* b_fact (- psiL psiR))
  (symbolic-simp '(+ (* (* b_fact (* c c)) (- BxL BxR)) (* 0.0 (- psiL psiR))))
   (symbolic-simp-rule '(+ (* (* b_fact (* c c)) (- BxL BxR)) (* 0.0 (- psiL psiR))))
    (symbolic-simp-rule '(* (* b_fact (* c c)) (- BxL BxR)))
    '(* b_fact (* (* c c) (- BxL BxR)))
    (symbolic-simp-rule '(* 0.0 (- psiL psiR)))
    0.0
   '(+ (* b_fact (* (* c c) (- BxL BxR))) 0.0)
  (symbolic-simp '(+ (* b_fact (* (* c c) (- BxL BxR))) 0.0))
   (symbolic-simp-rule '(+ (* b_fact (* (* c c) (- BxL BxR))) 0.0))
   '(+ 0.0 (* b_fact (* (* c c) (- BxL BxR))))
  (symbolic-simp '(+ 0.0 (* b_fact (* (* c c) (- BxL BxR)))))
   (symbolic-simp-rule '(+ 0.0 (* b_fact (* (* c c) (- BxL BxR)))))
   '(* b_fact (* (* c c) (- BxL BxR)))
  (symbolic-simp '(* b_fact (* (* c c) (- BxL BxR))))
   (symbolic-simp-rule '(* b_fact (* (* c c) (- BxL BxR))))
    (symbolic-simp-rule 'b_fact)
    'b_fact
    (symbolic-simp-rule '(* (* c c) (- BxL BxR)))
    '(* c (* c (- BxL BxR)))
   '(* b_fact (* c (* c (- BxL BxR))))
  (symbolic-simp '(* b_fact (* c (* c (- BxL BxR)))))
   (symbolic-simp-rule '(* b_fact (* c (* c (- BxL BxR)))))
    (symbolic-simp-rule 'b_fact)
    'b_fact
    (symbolic-simp-rule '(* c (* c (- BxL BxR))))
     (symbolic-simp-rule 'c)
     'c
     (symbolic-simp-rule '(* c (- BxL BxR)))
      (symbolic-simp-rule 'c)
      'c
      (symbolic-simp-rule '(- BxL BxR))
       (symbolic-simp-rule 'BxL)
       'BxL
       (symbolic-simp-rule 'BxR)
       'BxR
      '(- BxL BxR)
     '(* c (- BxL BxR))
    '(* c (* c (- BxL BxR)))
   '(* b_fact (* c (* c (- BxL BxR))))
  '(* b_fact (* c (* c (- BxL BxR))))
  (flux-deriv-replace '(* b_fact psi) 'Bx 'BxL)
   (flux-deriv-replace 'b_fact 'Bx 'BxL)
   'b_fact
   (flux-deriv-replace 'psi 'Bx 'BxL)
   'psi
  '(* b_fact psi)
  (flux-deriv-replace '(* b_fact psi) 'psi 'psiL)
   (flux-deriv-replace 'b_fact 'psi 'psiL)
   'b_fact
   (flux-deriv-replace 'psi 'psi 'psiL)
   'psiL
  '(* b_fact psiL)
  (flux-deriv-replace '(* b_fact psi) 'Bx 'BxR)
   (flux-deriv-replace 'b_fact 'Bx 'BxR)
   'b_fact
   (flux-deriv-replace 'psi 'Bx 'BxR)
   'psi
  '(* b_fact psi)
  (flux-deriv-replace '(* b_fact psi) 'psi 'psiR)
   (flux-deriv-replace 'b_fact 'psi 'psiR)
   'b_fact
   (flux-deriv-replace 'psi 'psi 'psiR)
   'psiR
  '(* b_fact psiR)
  (symbolic-simp '(- (* b_fact psiL) (* b_fact psiR)))
   (symbolic-simp-rule '(- (* b_fact psiL) (* b_fact psiR)))
   '(* b_fact (- psiL psiR))
  (symbolic-simp '(* b_fact (- psiL psiR)))
   (symbolic-simp-rule '(* b_fact (- psiL psiR)))
    (symbolic-simp-rule 'b_fact)
    'b_fact
    (symbolic-simp-rule '(- psiL psiR))
     (symbolic-simp-rule 'psiL)
     'psiL
     (symbolic-simp-rule 'psiR)
     'psiR
    '(- psiL psiR)
   '(* b_fact (- psiL psiR))
  '(* b_fact (- psiL psiR))
  (flux-deriv-replace '(* b_fact (* (* c c) Bx)) 'Bx 'BxL)
   (flux-deriv-replace 'b_fact 'Bx 'BxL)
   'b_fact
   (flux-deriv-replace '(* (* c c) Bx) 'Bx 'BxL)
    (flux-deriv-replace '(* c c) 'Bx 'BxL)
     (flux-deriv-replace 'c 'Bx 'BxL)
     'c
     (flux-deriv-replace 'c 'Bx 'BxL)
     'c
    '(* c c)
    (flux-deriv-replace 'Bx 'Bx 'BxL)
    'BxL
   '(* (* c c) BxL)
  '(* b_fact (* (* c c) BxL))
  (flux-deriv-replace '(* b_fact (* (* c c) BxL)) 'psi 'psiL)
   (flux-deriv-replace 'b_fact 'psi 'psiL)
   'b_fact
   (flux-deriv-replace '(* (* c c) BxL) 'psi 'psiL)
    (flux-deriv-replace '(* c c) 'psi 'psiL)
     (flux-deriv-replace 'c 'psi 'psiL)
     'c
     (flux-deriv-replace 'c 'psi 'psiL)
     'c
    '(* c c)
    (flux-deriv-replace 'BxL 'psi 'psiL)
    'BxL
   '(* (* c c) BxL)
  '(* b_fact (* (* c c) BxL))
  (flux-deriv-replace '(* b_fact (* (* c c) Bx)) 'Bx 'BxR)
   (flux-deriv-replace 'b_fact 'Bx 'BxR)
   'b_fact
   (flux-deriv-replace '(* (* c c) Bx) 'Bx 'BxR)
    (flux-deriv-replace '(* c c) 'Bx 'BxR)
     (flux-deriv-replace 'c 'Bx 'BxR)
     'c
     (flux-deriv-replace 'c 'Bx 'BxR)
     'c
    '(* c c)
    (flux-deriv-replace 'Bx 'Bx 'BxR)
    'BxR
   '(* (* c c) BxR)
  '(* b_fact (* (* c c) BxR))
  (flux-deriv-replace '(* b_fact (* (* c c) BxR)) 'psi 'psiR)
   (flux-deriv-replace 'b_fact 'psi 'psiR)
   'b_fact
   (flux-deriv-replace '(* (* c c) BxR) 'psi 'psiR)
    (flux-deriv-replace '(* c c) 'psi 'psiR)
     (flux-deriv-replace 'c 'psi 'psiR)
     'c
     (flux-deriv-replace 'c 'psi 'psiR)
     'c
    '(* c c)
    (flux-deriv-replace 'BxR 'psi 'psiR)
    'BxR
   '(* (* c c) BxR)
  '(* b_fact (* (* c c) BxR))
  (symbolic-simp '(- (* b_fact (* (* c c) BxL)) (* b_fact (* (* c c) BxR))))
   (symbolic-simp-rule '(- (* b_fact (* (* c c) BxL)) (* b_fact (* (* c c) BxR))))
   '(* b_fact (- (* (* c c) BxL) (* (* c c) BxR)))
  (symbolic-simp '(* b_fact (- (* (* c c) BxL) (* (* c c) BxR))))
   (symbolic-simp-rule '(* b_fact (- (* (* c c) BxL) (* (* c c) BxR))))
    (symbolic-simp-rule 'b_fact)
    'b_fact
    (symbolic-simp-rule '(- (* (* c c) BxL) (* (* c c) BxR)))
    '(* (* c c) (- BxL BxR))
   '(* b_fact (* (* c c) (- BxL BxR)))
  (symbolic-simp '(* b_fact (* (* c c) (- BxL BxR))))
   (symbolic-simp-rule '(* b_fact (* (* c c) (- BxL BxR))))
    (symbolic-simp-rule 'b_fact)
    'b_fact
    (symbolic-simp-rule '(* (* c c) (- BxL BxR)))
    '(* c (* c (- BxL BxR)))
   '(* b_fact (* c (* c (- BxL BxR))))
  (symbolic-simp '(* b_fact (* c (* c (- BxL BxR)))))
   (symbolic-simp-rule '(* b_fact (* c (* c (- BxL BxR)))))
    (symbolic-simp-rule 'b_fact)
    'b_fact
    (symbolic-simp-rule '(* c (* c (- BxL BxR))))
     (symbolic-simp-rule 'c)
     'c
     (symbolic-simp-rule '(* c (- BxL BxR)))
      (symbolic-simp-rule 'c)
      'c
      (symbolic-simp-rule '(- BxL BxR))
       (symbolic-simp-rule 'BxL)
       'BxL
       (symbolic-simp-rule 'BxR)
       'BxR
      '(- BxL BxR)
     '(* c (- BxL BxR))
    '(* c (* c (- BxL BxR)))
   '(* b_fact (* c (* c (- BxL BxR))))
  '(* b_fact (* c (* c (- BxL BxR))))
  (is-real 1.0 '(Bx psi) '((define c 1.0) (define e_fact 1.0) (define b_fact 1.0)))
  #t
  (is-real 1.0 '(Bx psi) '((define c 1.0) (define e_fact 1.0) (define b_fact 1.0)))
  #t
  (is-real 1.0 '(Bx psi) '((define c 1.0) (define e_fact 1.0) (define b_fact 1.0)))
  #t
  (is-real 0.0 '(Bx psi) '((define c 1.0) (define e_fact 1.0) (define b_fact 1.0)))
  #t
  (is-real '(cond ((< x 0.0) 0.5) (else -0.5)) '(Bx psi) '((define c 1.0) (define e_fact 1.0) (define b_fact 1.0)))
   (is-real 0.5 '(Bx psi) '((define c 1.0) (define e_fact 1.0) (define b_fact 1.0)))
   #t
  (is-real -0.5 '(Bx psi) '((define c 1.0) (define e_fact 1.0) (define b_fact 1.0)))
  #t
 #t
