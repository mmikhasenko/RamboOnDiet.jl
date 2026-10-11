# Edge cases

Above the total threshold, both deterministic maps use **open mass intervals and
closed angular intervals**:

- First `n - 2` inputs: strictly in `(0, 1)`.
- Remaining angular inputs: in `[0, 1]`.

An exact 0 or 1 in any mass coordinate raises `DomainError`. This policy is the
same for massive and massless daughters, even where a particular mass choice
admits a boundary limit. Mass endpoints collapse a sequential splitting; the
weight can diverge, or a massless cluster can lose its rest frame. The generic map
does not implement these boundary extensions. Angular endpoints remain valid,
including the fixed angles of aligned events. Two-body events have no mass inputs.

`generate_point`, `rand(generator)`, and `generate_aligned_point` redraw exact-zero
mass inputs internally, since Julia's floating-point `rand` samples `[0, 1)`.
The deterministic maps never redraw or clamp supplied inputs. Domain checks use
the working floating-point precision; a value that converts to an endpoint is
also rejected. No epsilon neighbourhood is excluded.

The existing whole-event threshold convention is separate: when the threshold
path is selected, finite inputs anywhere in the closed unit cube produce the
threshold configuration with zero weight. Nonfinite or out-of-cube inputs are
always rejected with `ArgumentError`.

For numerical integration above threshold, use interior nodes for the mass
coordinates. Excluding the boundary does not change the continuous phase-space
integral, but a quadrature rule that evaluates endpoints needs a separate boundary
extension or change of variables. Very small positive coordinates still require
appropriate numerical precision; the domain check is not an accuracy guarantee.

