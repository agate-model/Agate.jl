module Mortality

export linear_loss, quadratic_loss

linear_loss(P, l) = l * P

quadratic_loss(P, l) = l * P^2

end # module
