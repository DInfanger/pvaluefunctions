# representative error messages are stable

    Code
      cd_z(stderr = -1)
    Condition
      Error in `conf_dist()`:
      ! All values of `stderr` must be larger than 0.
      x You supplied -1.

---

    Code
      cd_z(together = NA)
    Condition
      Error in `conf_dist()`:
      ! `together` must be either `TRUE` or `FALSE`.
      x You supplied <logical> of length 1.

---

    Code
      cd_z(estimate = c(1, 2), stderr = c(1, 1), est_names = c("a", "a"))
    Condition
      Error in `conf_dist()`:
      ! Estimate names must be unique.
      x Duplicated name: "a".

