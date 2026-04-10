// traction from stress and normal
// --------------------------------------------------------------------------------------------
/** Calculate traction vector from stress for 2D plane strain elasticity.
 *
 * @param[in] stress Stress tensor.
 * @param[in] n Normal vector.
 * @param[out] traction Traction vector.
 */
static inline
void
traction(const pylith::fekernels::Tensor& stress,
         const PylithReal n[],
         PylithReal traction[]) {
    assert(traction);

    traction[0] = n[0]*stress.xx + n[1]*stress.xy;
    traction[1] = n[0]*stress.xy + n[1]*stress.yy;
} // traction
