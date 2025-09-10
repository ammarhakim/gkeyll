#lang racket

(require "../prover_core.rkt")

 (prove-flux-limiter-symmetry '#hash((limiter-expr . (max 0.0 (min 1.0 r))) (limiter-ratio . r) (name . "min_mod")))
  (symbolic-simp-positive '(/ (max 0.0 (min 1.0 r)) r) 'r)
   (symbolic-simp-positive-rule '(/ (max 0.0 (min 1.0 r)) r) 'r)
   '(max (/ (min 1.0 r) r) (/ 0.0 r))
  (symbolic-simp-positive '(max (/ (min 1.0 r) r) (/ 0.0 r)) 'r)
   (symbolic-simp-positive-rule '(max (/ (min 1.0 r) r) (/ 0.0 r)) 'r)
    (symbolic-simp-positive-rule '(/ (min 1.0 r) r) 'r)
    '(min (/ r r) (/ 1.0 r))
    (symbolic-simp-positive-rule '(/ 0.0 r) 'r)
     (symbolic-simp-positive-rule 0.0 'r)
     0.0
     (symbolic-simp-positive-rule 'r 'r)
     'r
    '(/ 0.0 r)
   '(max (min (/ r r) (/ 1.0 r)) (/ 0.0 r))
  (symbolic-simp-positive '(max (min (/ r r) (/ 1.0 r)) (/ 0.0 r)) 'r)
   (symbolic-simp-positive-rule '(max (min (/ r r) (/ 1.0 r)) (/ 0.0 r)) 'r)
    (symbolic-simp-positive-rule '(min (/ r r) (/ 1.0 r)) 'r)
     (symbolic-simp-positive-rule '(/ r r) 'r)
      (symbolic-simp-positive-rule 'r 'r)
      'r
      (symbolic-simp-positive-rule 'r 'r)
      'r
     '(/ r r)
     (symbolic-simp-positive-rule '(/ 1.0 r) 'r)
      (symbolic-simp-positive-rule 1.0 'r)
      1.0
      (symbolic-simp-positive-rule 'r 'r)
      'r
     '(/ 1.0 r)
    '(min (/ r r) (/ 1.0 r))
    (symbolic-simp-positive-rule '(/ 0.0 r) 'r)
     (symbolic-simp-positive-rule 0.0 'r)
     0.0
     (symbolic-simp-positive-rule 'r 'r)
     'r
    '(/ 0.0 r)
   '(max (min (/ r r) (/ 1.0 r)) (/ 0.0 r))
  '(max (min (/ r r) (/ 1.0 r)) (/ 0.0 r))
  (symbolic-simp '(max (min (/ r r) (/ 1.0 r)) (/ 0.0 r)))
   (symbolic-simp-rule '(max (min (/ r r) (/ 1.0 r)) (/ 0.0 r)))
   '(+ (* 0.5 (+ (min (/ r r) (/ 1.0 r)) (/ 0.0 r))) (* 0.5 (abs (- (min (/ r r) (/ 1.0 r)) (/ 0.0 r)))))
  (symbolic-simp '(+ (* 0.5 (+ (min (/ r r) (/ 1.0 r)) (/ 0.0 r))) (* 0.5 (abs (- (min (/ r r) (/ 1.0 r)) (/ 0.0 r))))))
   (symbolic-simp-rule '(+ (* 0.5 (+ (min (/ r r) (/ 1.0 r)) (/ 0.0 r))) (* 0.5 (abs (- (min (/ r r) (/ 1.0 r)) (/ 0.0 r))))))
    (symbolic-simp-rule '(* 0.5 (+ (min (/ r r) (/ 1.0 r)) (/ 0.0 r))))
    '(+ (* 0.5 (min (/ r r) (/ 1.0 r))) (* 0.5 (/ 0.0 r)))
    (symbolic-simp-rule '(* 0.5 (abs (- (min (/ r r) (/ 1.0 r)) (/ 0.0 r)))))
     (symbolic-simp-rule 0.5)
     0.5
     (symbolic-simp-rule '(abs (- (min (/ r r) (/ 1.0 r)) (/ 0.0 r))))
      (symbolic-simp-rule '(- (min (/ r r) (/ 1.0 r)) (/ 0.0 r)))
       (symbolic-simp-rule '(min (/ r r) (/ 1.0 r)))
       '(- (* 0.5 (+ (/ r r) (/ 1.0 r))) (* 0.5 (abs (- (/ r r) (/ 1.0 r)))))
       (symbolic-simp-rule '(/ 0.0 r))
       0.0
      '(- (- (* 0.5 (+ (/ r r) (/ 1.0 r))) (* 0.5 (abs (- (/ r r) (/ 1.0 r))))) 0.0)
     '(abs (- (- (* 0.5 (+ (/ r r) (/ 1.0 r))) (* 0.5 (abs (- (/ r r) (/ 1.0 r))))) 0.0))
    '(* 0.5 (abs (- (- (* 0.5 (+ (/ r r) (/ 1.0 r))) (* 0.5 (abs (- (/ r r) (/ 1.0 r))))) 0.0)))
   '(+ (+ (* 0.5 (min (/ r r) (/ 1.0 r))) (* 0.5 (/ 0.0 r))) (* 0.5 (abs (- (- (* 0.5 (+ (/ r r) (/ 1.0 r))) (* 0.5 (abs (- (/ r r) (/ 1.0 r))))) 0.0))))
  (symbolic-simp '(+ (+ (* 0.5 (min (/ r r) (/ 1.0 r))) (* 0.5 (/ 0.0 r))) (* 0.5 (abs (- (- (* 0.5 (+ (/ r r) (/ 1.0 r))) (* 0.5 (abs (- (/ r r) (/ 1.0 r))))) 0.0)))))
   (symbolic-simp-rule '(+ (+ (* 0.5 (min (/ r r) (/ 1.0 r))) (* 0.5 (/ 0.0 r))) (* 0.5 (abs (- (- (* 0.5 (+ (/ r r) (/ 1.0 r))) (* 0.5 (abs (- (/ r r) (/ 1.0 r))))) 0.0)))))
   '(+ (* 0.5 (min (/ r r) (/ 1.0 r))) (+ (* 0.5 (/ 0.0 r)) (* 0.5 (abs (- (- (* 0.5 (+ (/ r r) (/ 1.0 r))) (* 0.5 (abs (- (/ r r) (/ 1.0 r))))) 0.0)))))
  (symbolic-simp '(+ (* 0.5 (min (/ r r) (/ 1.0 r))) (+ (* 0.5 (/ 0.0 r)) (* 0.5 (abs (- (- (* 0.5 (+ (/ r r) (/ 1.0 r))) (* 0.5 (abs (- (/ r r) (/ 1.0 r))))) 0.0))))))
   (symbolic-simp-rule '(+ (* 0.5 (min (/ r r) (/ 1.0 r))) (+ (* 0.5 (/ 0.0 r)) (* 0.5 (abs (- (- (* 0.5 (+ (/ r r) (/ 1.0 r))) (* 0.5 (abs (- (/ r r) (/ 1.0 r))))) 0.0))))))
    (symbolic-simp-rule '(* 0.5 (min (/ r r) (/ 1.0 r))))
     (symbolic-simp-rule 0.5)
     0.5
     (symbolic-simp-rule '(min (/ r r) (/ 1.0 r)))
     '(- (* 0.5 (+ (/ r r) (/ 1.0 r))) (* 0.5 (abs (- (/ r r) (/ 1.0 r)))))
    '(* 0.5 (- (* 0.5 (+ (/ r r) (/ 1.0 r))) (* 0.5 (abs (- (/ r r) (/ 1.0 r))))))
    (symbolic-simp-rule '(+ (* 0.5 (/ 0.0 r)) (* 0.5 (abs (- (- (* 0.5 (+ (/ r r) (/ 1.0 r))) (* 0.5 (abs (- (/ r r) (/ 1.0 r))))) 0.0)))))
     (symbolic-simp-rule '(* 0.5 (/ 0.0 r)))
     '(/ 0.0 r)
     (symbolic-simp-rule '(* 0.5 (abs (- (- (* 0.5 (+ (/ r r) (/ 1.0 r))) (* 0.5 (abs (- (/ r r) (/ 1.0 r))))) 0.0))))
      (symbolic-simp-rule 0.5)
      0.5
      (symbolic-simp-rule '(abs (- (- (* 0.5 (+ (/ r r) (/ 1.0 r))) (* 0.5 (abs (- (/ r r) (/ 1.0 r))))) 0.0)))
       (symbolic-simp-rule '(- (- (* 0.5 (+ (/ r r) (/ 1.0 r))) (* 0.5 (abs (- (/ r r) (/ 1.0 r))))) 0.0))
       '(- (* 0.5 (+ (/ r r) (/ 1.0 r))) (* 0.5 (abs (- (/ r r) (/ 1.0 r)))))
      '(abs (- (* 0.5 (+ (/ r r) (/ 1.0 r))) (* 0.5 (abs (- (/ r r) (/ 1.0 r))))))
     '(* 0.5 (abs (- (* 0.5 (+ (/ r r) (/ 1.0 r))) (* 0.5 (abs (- (/ r r) (/ 1.0 r)))))))
    '(+ (/ 0.0 r) (* 0.5 (abs (- (* 0.5 (+ (/ r r) (/ 1.0 r))) (* 0.5 (abs (- (/ r r) (/ 1.0 r))))))))
   '(+ (* 0.5 (- (* 0.5 (+ (/ r r) (/ 1.0 r))) (* 0.5 (abs (- (/ r r) (/ 1.0 r)))))) (+ (/ 0.0 r) (* 0.5 (abs (- (* 0.5 (+ (/ r r) (/ 1.0 r))) (* 0.5 (abs (- (/ r r) (/ 1.0 r)))))))))
  (symbolic-simp '(+ (* 0.5 (- (* 0.5 (+ (/ r r) (/ 1.0 r))) (* 0.5 (abs (- (/ r r) (/ 1.0 r)))))) (+ (/ 0.0 r) (* 0.5 (abs (- (* 0.5 (+ (/ r r) (/ 1.0 r))) (* 0.5 (abs (- (/ r r) (/ 1.0 r))))))))))
   (symbolic-simp-rule '(+ (* 0.5 (- (* 0.5 (+ (/ r r) (/ 1.0 r))) (* 0.5 (abs (- (/ r r) (/ 1.0 r)))))) (+ (/ 0.0 r) (* 0.5 (abs (- (* 0.5 (+ (/ r r) (/ 1.0 r))) (* 0.5 (abs (- (/ r r) (/ 1.0 r))))))))))
    (symbolic-simp-rule '(* 0.5 (- (* 0.5 (+ (/ r r) (/ 1.0 r))) (* 0.5 (abs (- (/ r r) (/ 1.0 r)))))))
     (symbolic-simp-rule 0.5)
     0.5
     (symbolic-simp-rule '(- (* 0.5 (+ (/ r r) (/ 1.0 r))) (* 0.5 (abs (- (/ r r) (/ 1.0 r))))))
     '(* 0.5 (- (+ (/ r r) (/ 1.0 r)) (abs (- (/ r r) (/ 1.0 r)))))
    '(* 0.5 (* 0.5 (- (+ (/ r r) (/ 1.0 r)) (abs (- (/ r r) (/ 1.0 r))))))
    (symbolic-simp-rule '(+ (/ 0.0 r) (* 0.5 (abs (- (* 0.5 (+ (/ r r) (/ 1.0 r))) (* 0.5 (abs (- (/ r r) (/ 1.0 r)))))))))
     (symbolic-simp-rule '(/ 0.0 r))
     0.0
     (symbolic-simp-rule '(* 0.5 (abs (- (* 0.5 (+ (/ r r) (/ 1.0 r))) (* 0.5 (abs (- (/ r r) (/ 1.0 r))))))))
      (symbolic-simp-rule 0.5)
      0.5
      (symbolic-simp-rule '(abs (- (* 0.5 (+ (/ r r) (/ 1.0 r))) (* 0.5 (abs (- (/ r r) (/ 1.0 r)))))))
       (symbolic-simp-rule '(- (* 0.5 (+ (/ r r) (/ 1.0 r))) (* 0.5 (abs (- (/ r r) (/ 1.0 r))))))
       '(* 0.5 (- (+ (/ r r) (/ 1.0 r)) (abs (- (/ r r) (/ 1.0 r)))))
      '(abs (* 0.5 (- (+ (/ r r) (/ 1.0 r)) (abs (- (/ r r) (/ 1.0 r))))))
     '(* 0.5 (abs (* 0.5 (- (+ (/ r r) (/ 1.0 r)) (abs (- (/ r r) (/ 1.0 r)))))))
    '(+ 0.0 (* 0.5 (abs (* 0.5 (- (+ (/ r r) (/ 1.0 r)) (abs (- (/ r r) (/ 1.0 r))))))))
   '(+ (* 0.5 (* 0.5 (- (+ (/ r r) (/ 1.0 r)) (abs (- (/ r r) (/ 1.0 r)))))) (+ 0.0 (* 0.5 (abs (* 0.5 (- (+ (/ r r) (/ 1.0 r)) (abs (- (/ r r) (/ 1.0 r)))))))))
  (symbolic-simp '(+ (* 0.5 (* 0.5 (- (+ (/ r r) (/ 1.0 r)) (abs (- (/ r r) (/ 1.0 r)))))) (+ 0.0 (* 0.5 (abs (* 0.5 (- (+ (/ r r) (/ 1.0 r)) (abs (- (/ r r) (/ 1.0 r))))))))))
   (symbolic-simp-rule '(+ (* 0.5 (* 0.5 (- (+ (/ r r) (/ 1.0 r)) (abs (- (/ r r) (/ 1.0 r)))))) (+ 0.0 (* 0.5 (abs (* 0.5 (- (+ (/ r r) (/ 1.0 r)) (abs (- (/ r r) (/ 1.0 r))))))))))
    (symbolic-simp-rule '(* 0.5 (* 0.5 (- (+ (/ r r) (/ 1.0 r)) (abs (- (/ r r) (/ 1.0 r)))))))
    '(* 0.25 (- (+ (/ r r) (/ 1.0 r)) (abs (- (/ r r) (/ 1.0 r)))))
    (symbolic-simp-rule '(+ 0.0 (* 0.5 (abs (* 0.5 (- (+ (/ r r) (/ 1.0 r)) (abs (- (/ r r) (/ 1.0 r)))))))))
    '(* 0.5 (abs (* 0.5 (- (+ (/ r r) (/ 1.0 r)) (abs (- (/ r r) (/ 1.0 r)))))))
   '(+ (* 0.25 (- (+ (/ r r) (/ 1.0 r)) (abs (- (/ r r) (/ 1.0 r))))) (* 0.5 (abs (* 0.5 (- (+ (/ r r) (/ 1.0 r)) (abs (- (/ r r) (/ 1.0 r))))))))
  (symbolic-simp '(+ (* 0.25 (- (+ (/ r r) (/ 1.0 r)) (abs (- (/ r r) (/ 1.0 r))))) (* 0.5 (abs (* 0.5 (- (+ (/ r r) (/ 1.0 r)) (abs (- (/ r r) (/ 1.0 r)))))))))
   (symbolic-simp-rule '(+ (* 0.25 (- (+ (/ r r) (/ 1.0 r)) (abs (- (/ r r) (/ 1.0 r))))) (* 0.5 (abs (* 0.5 (- (+ (/ r r) (/ 1.0 r)) (abs (- (/ r r) (/ 1.0 r)))))))))
    (symbolic-simp-rule '(* 0.25 (- (+ (/ r r) (/ 1.0 r)) (abs (- (/ r r) (/ 1.0 r))))))
     (symbolic-simp-rule 0.25)
     0.25
     (symbolic-simp-rule '(- (+ (/ r r) (/ 1.0 r)) (abs (- (/ r r) (/ 1.0 r)))))
      (symbolic-simp-rule '(+ (/ r r) (/ 1.0 r)))
       (symbolic-simp-rule '(/ r r))
       1.0
       (symbolic-simp-rule '(/ 1.0 r))
        (symbolic-simp-rule 1.0)
        1.0
        (symbolic-simp-rule 'r)
        'r
       '(/ 1.0 r)
      '(+ 1.0 (/ 1.0 r))
      (symbolic-simp-rule '(abs (- (/ r r) (/ 1.0 r))))
       (symbolic-simp-rule '(- (/ r r) (/ 1.0 r)))
        (symbolic-simp-rule '(/ r r))
        1.0
        (symbolic-simp-rule '(/ 1.0 r))
         (symbolic-simp-rule 1.0)
         1.0
         (symbolic-simp-rule 'r)
         'r
        '(/ 1.0 r)
       '(- 1.0 (/ 1.0 r))
      '(abs (- 1.0 (/ 1.0 r)))
     '(- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r))))
    '(* 0.25 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r)))))
    (symbolic-simp-rule '(* 0.5 (abs (* 0.5 (- (+ (/ r r) (/ 1.0 r)) (abs (- (/ r r) (/ 1.0 r))))))))
     (symbolic-simp-rule 0.5)
     0.5
     (symbolic-simp-rule '(abs (* 0.5 (- (+ (/ r r) (/ 1.0 r)) (abs (- (/ r r) (/ 1.0 r)))))))
      (symbolic-simp-rule '(* 0.5 (- (+ (/ r r) (/ 1.0 r)) (abs (- (/ r r) (/ 1.0 r))))))
       (symbolic-simp-rule 0.5)
       0.5
       (symbolic-simp-rule '(- (+ (/ r r) (/ 1.0 r)) (abs (- (/ r r) (/ 1.0 r)))))
        (symbolic-simp-rule '(+ (/ r r) (/ 1.0 r)))
         (symbolic-simp-rule '(/ r r))
         1.0
         (symbolic-simp-rule '(/ 1.0 r))
          (symbolic-simp-rule 1.0)
          1.0
          (symbolic-simp-rule 'r)
          'r
         '(/ 1.0 r)
        '(+ 1.0 (/ 1.0 r))
        (symbolic-simp-rule '(abs (- (/ r r) (/ 1.0 r))))
         (symbolic-simp-rule '(- (/ r r) (/ 1.0 r)))
          (symbolic-simp-rule '(/ r r))
          1.0
          (symbolic-simp-rule '(/ 1.0 r))
        (symbolic-simp-rule 1.0)
        1.0
        (symbolic-simp-rule 'r)
        'r
          '(/ 1.0 r)
         '(- 1.0 (/ 1.0 r))
        '(abs (- 1.0 (/ 1.0 r)))
       '(- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r))))
      '(* 0.5 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r)))))
     '(abs (* 0.5 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r))))))
    '(* 0.5 (abs (* 0.5 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r)))))))
   '(+ (* 0.25 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r))))) (* 0.5 (abs (* 0.5 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r))))))))
  (symbolic-simp '(+ (* 0.25 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r))))) (* 0.5 (abs (* 0.5 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r)))))))))
   (symbolic-simp-rule '(+ (* 0.25 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r))))) (* 0.5 (abs (* 0.5 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r)))))))))
    (symbolic-simp-rule '(* 0.25 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r))))))
     (symbolic-simp-rule 0.25)
     0.25
     (symbolic-simp-rule '(- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r)))))
      (symbolic-simp-rule '(+ 1.0 (/ 1.0 r)))
       (symbolic-simp-rule 1.0)
       1.0
       (symbolic-simp-rule '(/ 1.0 r))
        (symbolic-simp-rule 1.0)
        1.0
        (symbolic-simp-rule 'r)
        'r
       '(/ 1.0 r)
      '(+ 1.0 (/ 1.0 r))
      (symbolic-simp-rule '(abs (- 1.0 (/ 1.0 r))))
       (symbolic-simp-rule '(- 1.0 (/ 1.0 r)))
        (symbolic-simp-rule 1.0)
        1.0
        (symbolic-simp-rule '(/ 1.0 r))
         (symbolic-simp-rule 1.0)
         1.0
         (symbolic-simp-rule 'r)
         'r
        '(/ 1.0 r)
       '(- 1.0 (/ 1.0 r))
      '(abs (- 1.0 (/ 1.0 r)))
     '(- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r))))
    '(* 0.25 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r)))))
    (symbolic-simp-rule '(* 0.5 (abs (* 0.5 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r))))))))
     (symbolic-simp-rule 0.5)
     0.5
     (symbolic-simp-rule '(abs (* 0.5 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r)))))))
      (symbolic-simp-rule '(* 0.5 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r))))))
       (symbolic-simp-rule 0.5)
       0.5
       (symbolic-simp-rule '(- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r)))))
        (symbolic-simp-rule '(+ 1.0 (/ 1.0 r)))
         (symbolic-simp-rule 1.0)
         1.0
         (symbolic-simp-rule '(/ 1.0 r))
          (symbolic-simp-rule 1.0)
          1.0
          (symbolic-simp-rule 'r)
          'r
         '(/ 1.0 r)
        '(+ 1.0 (/ 1.0 r))
        (symbolic-simp-rule '(abs (- 1.0 (/ 1.0 r))))
         (symbolic-simp-rule '(- 1.0 (/ 1.0 r)))
          (symbolic-simp-rule 1.0)
          1.0
          (symbolic-simp-rule '(/ 1.0 r))
        (symbolic-simp-rule 1.0)
        1.0
        (symbolic-simp-rule 'r)
        'r
          '(/ 1.0 r)
         '(- 1.0 (/ 1.0 r))
        '(abs (- 1.0 (/ 1.0 r)))
       '(- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r))))
      '(* 0.5 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r)))))
     '(abs (* 0.5 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r))))))
    '(* 0.5 (abs (* 0.5 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r)))))))
   '(+ (* 0.25 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r))))) (* 0.5 (abs (* 0.5 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r))))))))
  '(+ (* 0.25 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r))))) (* 0.5 (abs (* 0.5 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r))))))))
  (symbolic-simp-positive '(+ (* 0.25 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r))))) (* 0.5 (abs (* 0.5 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r)))))))) 'r)
   (symbolic-simp-positive-rule '(+ (* 0.25 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r))))) (* 0.5 (abs (* 0.5 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r)))))))) 'r)
    (symbolic-simp-positive-rule '(* 0.25 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r))))) 'r)
     (symbolic-simp-positive-rule 0.25 'r)
     0.25
     (symbolic-simp-positive-rule '(- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r)))) 'r)
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
      (symbolic-simp-positive-rule '(abs (- 1.0 (/ 1.0 r))) 'r)
      '(- 1.0 (/ 1.0 r))
     '(- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))
    '(* 0.25 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r))))
    (symbolic-simp-positive-rule '(* 0.5 (abs (* 0.5 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r))))))) 'r)
     (symbolic-simp-positive-rule 0.5 'r)
     0.5
     (symbolic-simp-positive-rule '(abs (* 0.5 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r)))))) 'r)
     '(* 0.5 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r)))))
    '(* 0.5 (* 0.5 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r))))))
   '(+ (* 0.25 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))) (* 0.5 (* 0.5 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r)))))))
  (symbolic-simp-positive '(+ (* 0.25 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))) (* 0.5 (* 0.5 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r))))))) 'r)
   (symbolic-simp-positive-rule '(+ (* 0.25 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))) (* 0.5 (* 0.5 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r))))))) 'r)
    (symbolic-simp-positive-rule '(* 0.25 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))) 'r)
     (symbolic-simp-positive-rule 0.25 'r)
     0.25
     (symbolic-simp-positive-rule '(- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r))) 'r)
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
      (symbolic-simp-positive-rule '(- 1.0 (/ 1.0 r)) 'r)
       (symbolic-simp-positive-rule 1.0 'r)
       1.0
       (symbolic-simp-positive-rule '(/ 1.0 r) 'r)
        (symbolic-simp-positive-rule 1.0 'r)
        1.0
        (symbolic-simp-positive-rule 'r 'r)
        'r
       '(/ 1.0 r)
      '(- 1.0 (/ 1.0 r))
     '(- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))
    '(* 0.25 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r))))
    (symbolic-simp-positive-rule '(* 0.5 (* 0.5 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r)))))) 'r)
     (symbolic-simp-positive-rule 0.5 'r)
     0.5
     (symbolic-simp-positive-rule '(* 0.5 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r))))) 'r)
      (symbolic-simp-positive-rule 0.5 'r)
      0.5
      (symbolic-simp-positive-rule '(- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r)))) 'r)
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
       (symbolic-simp-positive-rule '(abs (- 1.0 (/ 1.0 r))) 'r)
       '(- 1.0 (/ 1.0 r))
      '(- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))
     '(* 0.5 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r))))
    '(* 0.5 (* 0.5 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))))
   '(+ (* 0.25 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))) (* 0.5 (* 0.5 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r))))))
  (symbolic-simp-positive '(+ (* 0.25 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))) (* 0.5 (* 0.5 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))))) 'r)
   (symbolic-simp-positive-rule '(+ (* 0.25 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))) (* 0.5 (* 0.5 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))))) 'r)
    (symbolic-simp-positive-rule '(* 0.25 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))) 'r)
     (symbolic-simp-positive-rule 0.25 'r)
     0.25
     (symbolic-simp-positive-rule '(- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r))) 'r)
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
      (symbolic-simp-positive-rule '(- 1.0 (/ 1.0 r)) 'r)
       (symbolic-simp-positive-rule 1.0 'r)
       1.0
       (symbolic-simp-positive-rule '(/ 1.0 r) 'r)
        (symbolic-simp-positive-rule 1.0 'r)
        1.0
        (symbolic-simp-positive-rule 'r 'r)
        'r
       '(/ 1.0 r)
      '(- 1.0 (/ 1.0 r))
     '(- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))
    '(* 0.25 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r))))
    (symbolic-simp-positive-rule '(* 0.5 (* 0.5 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r))))) 'r)
     (symbolic-simp-positive-rule 0.5 'r)
     0.5
     (symbolic-simp-positive-rule '(* 0.5 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))) 'r)
      (symbolic-simp-positive-rule 0.5 'r)
      0.5
      (symbolic-simp-positive-rule '(- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r))) 'r)
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
       (symbolic-simp-positive-rule '(- 1.0 (/ 1.0 r)) 'r)
        (symbolic-simp-positive-rule 1.0 'r)
        1.0
        (symbolic-simp-positive-rule '(/ 1.0 r) 'r)
         (symbolic-simp-positive-rule 1.0 'r)
         1.0
         (symbolic-simp-positive-rule 'r 'r)
         'r
        '(/ 1.0 r)
       '(- 1.0 (/ 1.0 r))
      '(- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))
     '(* 0.5 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r))))
    '(* 0.5 (* 0.5 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))))
   '(+ (* 0.25 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))) (* 0.5 (* 0.5 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r))))))
  '(+ (* 0.25 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))) (* 0.5 (* 0.5 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r))))))
  (symbolic-simp '(+ (* 0.25 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))) (* 0.5 (* 0.5 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))))))
   (symbolic-simp-rule '(+ (* 0.25 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))) (* 0.5 (* 0.5 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))))))
    (symbolic-simp-rule '(* 0.25 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))))
     (symbolic-simp-rule 0.25)
     0.25
     (symbolic-simp-rule '(- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r))))
      (symbolic-simp-rule '(+ 1.0 (/ 1.0 r)))
       (symbolic-simp-rule 1.0)
       1.0
       (symbolic-simp-rule '(/ 1.0 r))
        (symbolic-simp-rule 1.0)
        1.0
        (symbolic-simp-rule 'r)
        'r
       '(/ 1.0 r)
      '(+ 1.0 (/ 1.0 r))
      (symbolic-simp-rule '(- 1.0 (/ 1.0 r)))
       (symbolic-simp-rule 1.0)
       1.0
       (symbolic-simp-rule '(/ 1.0 r))
        (symbolic-simp-rule 1.0)
        1.0
        (symbolic-simp-rule 'r)
        'r
       '(/ 1.0 r)
      '(- 1.0 (/ 1.0 r))
     '(- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))
    '(* 0.25 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r))))
    (symbolic-simp-rule '(* 0.5 (* 0.5 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r))))))
    '(* 0.25 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r))))
   '(+ (* 0.25 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))) (* 0.25 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))))
  (symbolic-simp '(+ (* 0.25 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))) (* 0.25 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r))))))
   (symbolic-simp-rule '(+ (* 0.25 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))) (* 0.25 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r))))))
   '(* (+ 0.25 0.25) (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r))))
  (symbolic-simp '(* (+ 0.25 0.25) (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))))
   (symbolic-simp-rule '(* (+ 0.25 0.25) (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))))
    (symbolic-simp-rule '(+ 0.25 0.25))
    0.5
    (symbolic-simp-rule '(- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r))))
     (symbolic-simp-rule '(+ 1.0 (/ 1.0 r)))
      (symbolic-simp-rule 1.0)
      1.0
      (symbolic-simp-rule '(/ 1.0 r))
       (symbolic-simp-rule 1.0)
       1.0
       (symbolic-simp-rule 'r)
       'r
      '(/ 1.0 r)
     '(+ 1.0 (/ 1.0 r))
     (symbolic-simp-rule '(- 1.0 (/ 1.0 r)))
      (symbolic-simp-rule 1.0)
      1.0
      (symbolic-simp-rule '(/ 1.0 r))
       (symbolic-simp-rule 1.0)
       1.0
       (symbolic-simp-rule 'r)
       'r
      '(/ 1.0 r)
     '(- 1.0 (/ 1.0 r))
    '(- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))
   '(* 0.5 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r))))
  (symbolic-simp '(* 0.5 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))))
   (symbolic-simp-rule '(* 0.5 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))))
    (symbolic-simp-rule 0.5)
    0.5
    (symbolic-simp-rule '(- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r))))
     (symbolic-simp-rule '(+ 1.0 (/ 1.0 r)))
      (symbolic-simp-rule 1.0)
      1.0
      (symbolic-simp-rule '(/ 1.0 r))
       (symbolic-simp-rule 1.0)
       1.0
       (symbolic-simp-rule 'r)
       'r
      '(/ 1.0 r)
     '(+ 1.0 (/ 1.0 r))
     (symbolic-simp-rule '(- 1.0 (/ 1.0 r)))
      (symbolic-simp-rule 1.0)
      1.0
      (symbolic-simp-rule '(/ 1.0 r))
       (symbolic-simp-rule 1.0)
       1.0
       (symbolic-simp-rule 'r)
       'r
      '(/ 1.0 r)
     '(- 1.0 (/ 1.0 r))
    '(- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))
   '(* 0.5 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r))))
  '(* 0.5 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r))))
  (variable-transform '(max 0.0 (min 1.0 r)) 'r '(/ 1.0 r))
   (variable-transform 'max 'r '(/ 1.0 r))
   'max
   (variable-transform 0.0 'r '(/ 1.0 r))
   0.0
   (variable-transform '(min 1.0 r) 'r '(/ 1.0 r))
    (variable-transform 'min 'r '(/ 1.0 r))
    'min
    (variable-transform 1.0 'r '(/ 1.0 r))
    1.0
    (variable-transform 'r 'r '(/ 1.0 r))
    '(/ 1.0 r)
   '(min 1.0 (/ 1.0 r))
  '(max 0.0 (min 1.0 (/ 1.0 r)))
  (symbolic-simp-positive '(max 0.0 (min 1.0 (/ 1.0 r))) 'r)
   (symbolic-simp-positive-rule '(max 0.0 (min 1.0 (/ 1.0 r))) 'r)
    (symbolic-simp-positive-rule 0.0 'r)
    0.0
    (symbolic-simp-positive-rule '(min 1.0 (/ 1.0 r)) 'r)
     (symbolic-simp-positive-rule 1.0 'r)
     1.0
     (symbolic-simp-positive-rule '(/ 1.0 r) 'r)
      (symbolic-simp-positive-rule 1.0 'r)
      1.0
      (symbolic-simp-positive-rule 'r 'r)
      'r
     '(/ 1.0 r)
    '(min 1.0 (/ 1.0 r))
   '(max 0.0 (min 1.0 (/ 1.0 r)))
  '(max 0.0 (min 1.0 (/ 1.0 r)))
  (symbolic-simp '(max 0.0 (min 1.0 (/ 1.0 r))))
   (symbolic-simp-rule '(max 0.0 (min 1.0 (/ 1.0 r))))
   '(+ (* 0.5 (+ 0.0 (min 1.0 (/ 1.0 r)))) (* 0.5 (abs (- 0.0 (min 1.0 (/ 1.0 r))))))
  (symbolic-simp '(+ (* 0.5 (+ 0.0 (min 1.0 (/ 1.0 r)))) (* 0.5 (abs (- 0.0 (min 1.0 (/ 1.0 r)))))))
   (symbolic-simp-rule '(+ (* 0.5 (+ 0.0 (min 1.0 (/ 1.0 r)))) (* 0.5 (abs (- 0.0 (min 1.0 (/ 1.0 r)))))))
    (symbolic-simp-rule '(* 0.5 (+ 0.0 (min 1.0 (/ 1.0 r)))))
    '(+ (* 0.5 0.0) (* 0.5 (min 1.0 (/ 1.0 r))))
    (symbolic-simp-rule '(* 0.5 (abs (- 0.0 (min 1.0 (/ 1.0 r))))))
     (symbolic-simp-rule 0.5)
     0.5
     (symbolic-simp-rule '(abs (- 0.0 (min 1.0 (/ 1.0 r)))))
      (symbolic-simp-rule '(- 0.0 (min 1.0 (/ 1.0 r))))
      '(* -1.0 (min 1.0 (/ 1.0 r)))
     '(abs (* -1.0 (min 1.0 (/ 1.0 r))))
    '(* 0.5 (abs (* -1.0 (min 1.0 (/ 1.0 r)))))
   '(+ (+ (* 0.5 0.0) (* 0.5 (min 1.0 (/ 1.0 r)))) (* 0.5 (abs (* -1.0 (min 1.0 (/ 1.0 r))))))
  (symbolic-simp '(+ (+ (* 0.5 0.0) (* 0.5 (min 1.0 (/ 1.0 r)))) (* 0.5 (abs (* -1.0 (min 1.0 (/ 1.0 r)))))))
   (symbolic-simp-rule '(+ (+ (* 0.5 0.0) (* 0.5 (min 1.0 (/ 1.0 r)))) (* 0.5 (abs (* -1.0 (min 1.0 (/ 1.0 r)))))))
   '(+ (* 0.5 0.0) (+ (* 0.5 (min 1.0 (/ 1.0 r))) (* 0.5 (abs (* -1.0 (min 1.0 (/ 1.0 r)))))))
  (symbolic-simp '(+ (* 0.5 0.0) (+ (* 0.5 (min 1.0 (/ 1.0 r))) (* 0.5 (abs (* -1.0 (min 1.0 (/ 1.0 r))))))))
   (symbolic-simp-rule '(+ (* 0.5 0.0) (+ (* 0.5 (min 1.0 (/ 1.0 r))) (* 0.5 (abs (* -1.0 (min 1.0 (/ 1.0 r))))))))
    (symbolic-simp-rule '(* 0.5 0.0))
    0.0
    (symbolic-simp-rule '(+ (* 0.5 (min 1.0 (/ 1.0 r))) (* 0.5 (abs (* -1.0 (min 1.0 (/ 1.0 r)))))))
     (symbolic-simp-rule '(* 0.5 (min 1.0 (/ 1.0 r))))
      (symbolic-simp-rule 0.5)
      0.5
      (symbolic-simp-rule '(min 1.0 (/ 1.0 r)))
      '(- (* 0.5 (+ 1.0 (/ 1.0 r))) (* 0.5 (abs (- 1.0 (/ 1.0 r)))))
     '(* 0.5 (- (* 0.5 (+ 1.0 (/ 1.0 r))) (* 0.5 (abs (- 1.0 (/ 1.0 r))))))
     (symbolic-simp-rule '(* 0.5 (abs (* -1.0 (min 1.0 (/ 1.0 r))))))
      (symbolic-simp-rule 0.5)
      0.5
      (symbolic-simp-rule '(abs (* -1.0 (min 1.0 (/ 1.0 r)))))
      '(abs (min 1.0 (/ 1.0 r)))
     '(* 0.5 (abs (min 1.0 (/ 1.0 r))))
    '(+ (* 0.5 (- (* 0.5 (+ 1.0 (/ 1.0 r))) (* 0.5 (abs (- 1.0 (/ 1.0 r)))))) (* 0.5 (abs (min 1.0 (/ 1.0 r)))))
   '(+ 0.0 (+ (* 0.5 (- (* 0.5 (+ 1.0 (/ 1.0 r))) (* 0.5 (abs (- 1.0 (/ 1.0 r)))))) (* 0.5 (abs (min 1.0 (/ 1.0 r))))))
  (symbolic-simp '(+ 0.0 (+ (* 0.5 (- (* 0.5 (+ 1.0 (/ 1.0 r))) (* 0.5 (abs (- 1.0 (/ 1.0 r)))))) (* 0.5 (abs (min 1.0 (/ 1.0 r)))))))
   (symbolic-simp-rule '(+ 0.0 (+ (* 0.5 (- (* 0.5 (+ 1.0 (/ 1.0 r))) (* 0.5 (abs (- 1.0 (/ 1.0 r)))))) (* 0.5 (abs (min 1.0 (/ 1.0 r)))))))
   '(+ (* 0.5 (- (* 0.5 (+ 1.0 (/ 1.0 r))) (* 0.5 (abs (- 1.0 (/ 1.0 r)))))) (* 0.5 (abs (min 1.0 (/ 1.0 r)))))
  (symbolic-simp '(+ (* 0.5 (- (* 0.5 (+ 1.0 (/ 1.0 r))) (* 0.5 (abs (- 1.0 (/ 1.0 r)))))) (* 0.5 (abs (min 1.0 (/ 1.0 r))))))
   (symbolic-simp-rule '(+ (* 0.5 (- (* 0.5 (+ 1.0 (/ 1.0 r))) (* 0.5 (abs (- 1.0 (/ 1.0 r)))))) (* 0.5 (abs (min 1.0 (/ 1.0 r))))))
    (symbolic-simp-rule '(* 0.5 (- (* 0.5 (+ 1.0 (/ 1.0 r))) (* 0.5 (abs (- 1.0 (/ 1.0 r)))))))
     (symbolic-simp-rule 0.5)
     0.5
     (symbolic-simp-rule '(- (* 0.5 (+ 1.0 (/ 1.0 r))) (* 0.5 (abs (- 1.0 (/ 1.0 r))))))
     '(* 0.5 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r)))))
    '(* 0.5 (* 0.5 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r))))))
    (symbolic-simp-rule '(* 0.5 (abs (min 1.0 (/ 1.0 r)))))
     (symbolic-simp-rule 0.5)
     0.5
     (symbolic-simp-rule '(abs (min 1.0 (/ 1.0 r))))
      (symbolic-simp-rule '(min 1.0 (/ 1.0 r)))
      '(- (* 0.5 (+ 1.0 (/ 1.0 r))) (* 0.5 (abs (- 1.0 (/ 1.0 r)))))
     '(abs (- (* 0.5 (+ 1.0 (/ 1.0 r))) (* 0.5 (abs (- 1.0 (/ 1.0 r))))))
    '(* 0.5 (abs (- (* 0.5 (+ 1.0 (/ 1.0 r))) (* 0.5 (abs (- 1.0 (/ 1.0 r)))))))
   '(+ (* 0.5 (* 0.5 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r)))))) (* 0.5 (abs (- (* 0.5 (+ 1.0 (/ 1.0 r))) (* 0.5 (abs (- 1.0 (/ 1.0 r))))))))
  (symbolic-simp '(+ (* 0.5 (* 0.5 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r)))))) (* 0.5 (abs (- (* 0.5 (+ 1.0 (/ 1.0 r))) (* 0.5 (abs (- 1.0 (/ 1.0 r)))))))))
   (symbolic-simp-rule '(+ (* 0.5 (* 0.5 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r)))))) (* 0.5 (abs (- (* 0.5 (+ 1.0 (/ 1.0 r))) (* 0.5 (abs (- 1.0 (/ 1.0 r)))))))))
    (symbolic-simp-rule '(* 0.5 (* 0.5 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r)))))))
    '(* 0.25 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r)))))
    (symbolic-simp-rule '(* 0.5 (abs (- (* 0.5 (+ 1.0 (/ 1.0 r))) (* 0.5 (abs (- 1.0 (/ 1.0 r))))))))
     (symbolic-simp-rule 0.5)
     0.5
     (symbolic-simp-rule '(abs (- (* 0.5 (+ 1.0 (/ 1.0 r))) (* 0.5 (abs (- 1.0 (/ 1.0 r)))))))
      (symbolic-simp-rule '(- (* 0.5 (+ 1.0 (/ 1.0 r))) (* 0.5 (abs (- 1.0 (/ 1.0 r))))))
      '(* 0.5 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r)))))
     '(abs (* 0.5 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r))))))
    '(* 0.5 (abs (* 0.5 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r)))))))
   '(+ (* 0.25 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r))))) (* 0.5 (abs (* 0.5 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r))))))))
  (symbolic-simp '(+ (* 0.25 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r))))) (* 0.5 (abs (* 0.5 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r)))))))))
   (symbolic-simp-rule '(+ (* 0.25 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r))))) (* 0.5 (abs (* 0.5 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r)))))))))
    (symbolic-simp-rule '(* 0.25 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r))))))
     (symbolic-simp-rule 0.25)
     0.25
     (symbolic-simp-rule '(- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r)))))
      (symbolic-simp-rule '(+ 1.0 (/ 1.0 r)))
       (symbolic-simp-rule 1.0)
       1.0
       (symbolic-simp-rule '(/ 1.0 r))
        (symbolic-simp-rule 1.0)
        1.0
        (symbolic-simp-rule 'r)
        'r
       '(/ 1.0 r)
      '(+ 1.0 (/ 1.0 r))
      (symbolic-simp-rule '(abs (- 1.0 (/ 1.0 r))))
       (symbolic-simp-rule '(- 1.0 (/ 1.0 r)))
        (symbolic-simp-rule 1.0)
        1.0
        (symbolic-simp-rule '(/ 1.0 r))
         (symbolic-simp-rule 1.0)
         1.0
         (symbolic-simp-rule 'r)
         'r
        '(/ 1.0 r)
       '(- 1.0 (/ 1.0 r))
      '(abs (- 1.0 (/ 1.0 r)))
     '(- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r))))
    '(* 0.25 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r)))))
    (symbolic-simp-rule '(* 0.5 (abs (* 0.5 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r))))))))
     (symbolic-simp-rule 0.5)
     0.5
     (symbolic-simp-rule '(abs (* 0.5 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r)))))))
      (symbolic-simp-rule '(* 0.5 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r))))))
       (symbolic-simp-rule 0.5)
       0.5
       (symbolic-simp-rule '(- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r)))))
        (symbolic-simp-rule '(+ 1.0 (/ 1.0 r)))
         (symbolic-simp-rule 1.0)
         1.0
         (symbolic-simp-rule '(/ 1.0 r))
          (symbolic-simp-rule 1.0)
          1.0
          (symbolic-simp-rule 'r)
          'r
         '(/ 1.0 r)
        '(+ 1.0 (/ 1.0 r))
        (symbolic-simp-rule '(abs (- 1.0 (/ 1.0 r))))
         (symbolic-simp-rule '(- 1.0 (/ 1.0 r)))
          (symbolic-simp-rule 1.0)
          1.0
          (symbolic-simp-rule '(/ 1.0 r))
        (symbolic-simp-rule 1.0)
        1.0
        (symbolic-simp-rule 'r)
        'r
          '(/ 1.0 r)
         '(- 1.0 (/ 1.0 r))
        '(abs (- 1.0 (/ 1.0 r)))
       '(- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r))))
      '(* 0.5 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r)))))
     '(abs (* 0.5 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r))))))
    '(* 0.5 (abs (* 0.5 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r)))))))
   '(+ (* 0.25 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r))))) (* 0.5 (abs (* 0.5 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r))))))))
  '(+ (* 0.25 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r))))) (* 0.5 (abs (* 0.5 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r))))))))
  (symbolic-simp-positive '(+ (* 0.25 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r))))) (* 0.5 (abs (* 0.5 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r)))))))) 'r)
   (symbolic-simp-positive-rule '(+ (* 0.25 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r))))) (* 0.5 (abs (* 0.5 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r)))))))) 'r)
    (symbolic-simp-positive-rule '(* 0.25 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r))))) 'r)
     (symbolic-simp-positive-rule 0.25 'r)
     0.25
     (symbolic-simp-positive-rule '(- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r)))) 'r)
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
      (symbolic-simp-positive-rule '(abs (- 1.0 (/ 1.0 r))) 'r)
      '(- 1.0 (/ 1.0 r))
     '(- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))
    '(* 0.25 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r))))
    (symbolic-simp-positive-rule '(* 0.5 (abs (* 0.5 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r))))))) 'r)
     (symbolic-simp-positive-rule 0.5 'r)
     0.5
     (symbolic-simp-positive-rule '(abs (* 0.5 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r)))))) 'r)
     '(* 0.5 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r)))))
    '(* 0.5 (* 0.5 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r))))))
   '(+ (* 0.25 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))) (* 0.5 (* 0.5 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r)))))))
  (symbolic-simp-positive '(+ (* 0.25 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))) (* 0.5 (* 0.5 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r))))))) 'r)
   (symbolic-simp-positive-rule '(+ (* 0.25 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))) (* 0.5 (* 0.5 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r))))))) 'r)
    (symbolic-simp-positive-rule '(* 0.25 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))) 'r)
     (symbolic-simp-positive-rule 0.25 'r)
     0.25
     (symbolic-simp-positive-rule '(- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r))) 'r)
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
      (symbolic-simp-positive-rule '(- 1.0 (/ 1.0 r)) 'r)
       (symbolic-simp-positive-rule 1.0 'r)
       1.0
       (symbolic-simp-positive-rule '(/ 1.0 r) 'r)
        (symbolic-simp-positive-rule 1.0 'r)
        1.0
        (symbolic-simp-positive-rule 'r 'r)
        'r
       '(/ 1.0 r)
      '(- 1.0 (/ 1.0 r))
     '(- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))
    '(* 0.25 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r))))
    (symbolic-simp-positive-rule '(* 0.5 (* 0.5 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r)))))) 'r)
     (symbolic-simp-positive-rule 0.5 'r)
     0.5
     (symbolic-simp-positive-rule '(* 0.5 (- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r))))) 'r)
      (symbolic-simp-positive-rule 0.5 'r)
      0.5
      (symbolic-simp-positive-rule '(- (+ 1.0 (/ 1.0 r)) (abs (- 1.0 (/ 1.0 r)))) 'r)
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
       (symbolic-simp-positive-rule '(abs (- 1.0 (/ 1.0 r))) 'r)
       '(- 1.0 (/ 1.0 r))
      '(- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))
     '(* 0.5 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r))))
    '(* 0.5 (* 0.5 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))))
   '(+ (* 0.25 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))) (* 0.5 (* 0.5 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r))))))
  (symbolic-simp-positive '(+ (* 0.25 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))) (* 0.5 (* 0.5 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))))) 'r)
   (symbolic-simp-positive-rule '(+ (* 0.25 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))) (* 0.5 (* 0.5 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))))) 'r)
    (symbolic-simp-positive-rule '(* 0.25 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))) 'r)
     (symbolic-simp-positive-rule 0.25 'r)
     0.25
     (symbolic-simp-positive-rule '(- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r))) 'r)
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
      (symbolic-simp-positive-rule '(- 1.0 (/ 1.0 r)) 'r)
       (symbolic-simp-positive-rule 1.0 'r)
       1.0
       (symbolic-simp-positive-rule '(/ 1.0 r) 'r)
        (symbolic-simp-positive-rule 1.0 'r)
        1.0
        (symbolic-simp-positive-rule 'r 'r)
        'r
       '(/ 1.0 r)
      '(- 1.0 (/ 1.0 r))
     '(- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))
    '(* 0.25 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r))))
    (symbolic-simp-positive-rule '(* 0.5 (* 0.5 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r))))) 'r)
     (symbolic-simp-positive-rule 0.5 'r)
     0.5
     (symbolic-simp-positive-rule '(* 0.5 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))) 'r)
      (symbolic-simp-positive-rule 0.5 'r)
      0.5
      (symbolic-simp-positive-rule '(- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r))) 'r)
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
       (symbolic-simp-positive-rule '(- 1.0 (/ 1.0 r)) 'r)
        (symbolic-simp-positive-rule 1.0 'r)
        1.0
        (symbolic-simp-positive-rule '(/ 1.0 r) 'r)
         (symbolic-simp-positive-rule 1.0 'r)
         1.0
         (symbolic-simp-positive-rule 'r 'r)
         'r
        '(/ 1.0 r)
       '(- 1.0 (/ 1.0 r))
      '(- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))
     '(* 0.5 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r))))
    '(* 0.5 (* 0.5 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))))
   '(+ (* 0.25 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))) (* 0.5 (* 0.5 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r))))))
  '(+ (* 0.25 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))) (* 0.5 (* 0.5 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r))))))
  (symbolic-simp '(+ (* 0.25 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))) (* 0.5 (* 0.5 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))))))
   (symbolic-simp-rule '(+ (* 0.25 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))) (* 0.5 (* 0.5 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))))))
    (symbolic-simp-rule '(* 0.25 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))))
     (symbolic-simp-rule 0.25)
     0.25
     (symbolic-simp-rule '(- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r))))
      (symbolic-simp-rule '(+ 1.0 (/ 1.0 r)))
       (symbolic-simp-rule 1.0)
       1.0
       (symbolic-simp-rule '(/ 1.0 r))
        (symbolic-simp-rule 1.0)
        1.0
        (symbolic-simp-rule 'r)
        'r
       '(/ 1.0 r)
      '(+ 1.0 (/ 1.0 r))
      (symbolic-simp-rule '(- 1.0 (/ 1.0 r)))
       (symbolic-simp-rule 1.0)
       1.0
       (symbolic-simp-rule '(/ 1.0 r))
        (symbolic-simp-rule 1.0)
        1.0
        (symbolic-simp-rule 'r)
        'r
       '(/ 1.0 r)
      '(- 1.0 (/ 1.0 r))
     '(- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))
    '(* 0.25 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r))))
    (symbolic-simp-rule '(* 0.5 (* 0.5 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r))))))
    '(* 0.25 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r))))
   '(+ (* 0.25 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))) (* 0.25 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))))
  (symbolic-simp '(+ (* 0.25 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))) (* 0.25 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r))))))
   (symbolic-simp-rule '(+ (* 0.25 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))) (* 0.25 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r))))))
   '(* (+ 0.25 0.25) (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r))))
  (symbolic-simp '(* (+ 0.25 0.25) (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))))
   (symbolic-simp-rule '(* (+ 0.25 0.25) (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))))
    (symbolic-simp-rule '(+ 0.25 0.25))
    0.5
    (symbolic-simp-rule '(- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r))))
     (symbolic-simp-rule '(+ 1.0 (/ 1.0 r)))
      (symbolic-simp-rule 1.0)
      1.0
      (symbolic-simp-rule '(/ 1.0 r))
       (symbolic-simp-rule 1.0)
       1.0
       (symbolic-simp-rule 'r)
       'r
      '(/ 1.0 r)
     '(+ 1.0 (/ 1.0 r))
     (symbolic-simp-rule '(- 1.0 (/ 1.0 r)))
      (symbolic-simp-rule 1.0)
      1.0
      (symbolic-simp-rule '(/ 1.0 r))
       (symbolic-simp-rule 1.0)
       1.0
       (symbolic-simp-rule 'r)
       'r
      '(/ 1.0 r)
     '(- 1.0 (/ 1.0 r))
    '(- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))
   '(* 0.5 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r))))
  (symbolic-simp '(* 0.5 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))))
   (symbolic-simp-rule '(* 0.5 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))))
    (symbolic-simp-rule 0.5)
    0.5
    (symbolic-simp-rule '(- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r))))
     (symbolic-simp-rule '(+ 1.0 (/ 1.0 r)))
      (symbolic-simp-rule 1.0)
      1.0
      (symbolic-simp-rule '(/ 1.0 r))
       (symbolic-simp-rule 1.0)
       1.0
       (symbolic-simp-rule 'r)
       'r
      '(/ 1.0 r)
     '(+ 1.0 (/ 1.0 r))
     (symbolic-simp-rule '(- 1.0 (/ 1.0 r)))
      (symbolic-simp-rule 1.0)
      1.0
      (symbolic-simp-rule '(/ 1.0 r))
       (symbolic-simp-rule 1.0)
       1.0
       (symbolic-simp-rule 'r)
       'r
      '(/ 1.0 r)
     '(- 1.0 (/ 1.0 r))
    '(- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r)))
   '(* 0.5 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r))))
  '(* 0.5 (- (+ 1.0 (/ 1.0 r)) (- 1.0 (/ 1.0 r))))
 #t
