#pragma once

namespace xmath
{

    //------------------------------------------------------------------------------
    // radian3
    //------------------------------------------------------------------------------
    //
    // Constructs a radian3 with all angles set to the same value.
    //
    // Parameters:
    //  Angles - The angle value to set for pitch, yaw, and roll.
    //
    constexpr radian3::radian3(const radian Angles) noexcept : m_Pitch(Angles), m_Yaw(Angles), m_Roll(Angles)
    {
      //  assert(xmath::isValid(m_Pitch.m_Value) && "Invalid angle value in constructor");
      //  assert(xmath::isValid(m_Yaw.m_Value)   && "Invalid angle value in constructor");
      //  assert(xmath::isValid(m_Roll.m_Value)  && "Invalid angle value in constructor");
    }

    //------------------------------------------------------------------------------
    // radian3
    //------------------------------------------------------------------------------
    //
    // Constructs a radian3 with specified pitch, yaw, and roll.
    //
    // Parameters:
    //  Pitch - The pitch angle.
    //  Yaw - The yaw angle.
    //  Roll - The roll angle.
    //
    constexpr radian3::radian3(radian Pitch, radian Yaw, radian Roll) noexcept : m_Pitch(Pitch), m_Yaw(Yaw), m_Roll(Roll)
    {
     //   assert(xmath::isValid(m_Pitch.m_Value) && "Invalid angle value in constructor");
     //   assert(xmath::isValid(m_Yaw.m_Value) && "Invalid angle value in constructor");
     //   assert(xmath::isValid(m_Roll.m_Value) && "Invalid angle value in constructor");
    }

    //------------------------------------------------------------------------------
    // fromZero
    //------------------------------------------------------------------------------
    //
    // Sets all angles to zero.
    //
    // Returns:
    //  Void.
    //
    inline void radian3::Zero(void) noexcept
    {
        m_Pitch = radian{ 0.0f };
        m_Yaw = radian{ 0.0f };
        m_Roll = radian{ 0.0f };
    }

    //------------------------------------------------------------------------------

    consteval radian3 radian3::fromZero(void) noexcept
    {
        return radian3{ radian{0}, radian{ 0 }, radian{0} };
    }

    //------------------------------------------------------------------------------
    // ModAngle
    //------------------------------------------------------------------------------
    //
    // Normalizes each angle to [0, 2pi).
    //
    // Returns:
    //  Reference to self.
    //
    inline radian3& radian3::ModAngle(void) noexcept
    {
        assert(isValid() && "Invalid angles before ModAngle");
        m_Pitch = xmath::ModAngle(m_Pitch);
        m_Yaw   = xmath::ModAngle(m_Yaw);
        m_Roll  = xmath::ModAngle(m_Roll);
        return *this;
    }

    //------------------------------------------------------------------------------
    // ModAngleCopy
    //------------------------------------------------------------------------------
    //
    // Returns a copy normalized to [0, 2pi).
    //
    // Returns:
    //  Normalized copy.
    //
    inline radian3 radian3::ModAngleCopy(void) const noexcept
    {
        assert(isValid() && "Invalid angles in ModAngleCopy");
        radian3 copy = *this;
        return copy.ModAngle();
    }

    //------------------------------------------------------------------------------
    // ModAngle2
    //------------------------------------------------------------------------------
    //
    // Normalizes each angle to [-pi, pi).
    //
    // Returns:
    //  Reference to self.
    //
    inline radian3& radian3::ModAngle2(void) noexcept
    {
        assert(isValid() && "Invalid angles before ModAngle2");
        m_Pitch = xmath::ModAngle2(m_Pitch);
        m_Yaw   = xmath::ModAngle2(m_Yaw);
        m_Roll  = xmath::ModAngle2(m_Roll);
        return *this;
    }

    //------------------------------------------------------------------------------
    // ModAngle2Copy
    //------------------------------------------------------------------------------
    //
    // Returns a copy normalized to [-pi, pi).
    //
    // Returns:
    //  Normalized copy.
    //
    inline radian3 radian3::ModAngle2Copy(void) const noexcept
    {
        assert(isValid() && "Invalid angles in ModAngle2Copy");
        radian3 copy = *this;
        return copy.ModAngle2();
    }

    //------------------------------------------------------------------------------
    // MinAngleDiff
    //------------------------------------------------------------------------------
    //
    // Computes the minimal angle difference for each component.
    //
    // Parameters:
    //  R - The other radian3 to compare with.
    //
    // Returns:
    //  radian3 with the minimal differences per component.
    //
    inline radian3 radian3::MinAngleDiff(const radian3& R) const noexcept
    {
        assert(isValid() && R.isValid() && "Invalid angles in MinAngleDiff");
        return radian3{
            xmath::MinAngleDiff(m_Pitch, R.m_Pitch),
            xmath::MinAngleDiff(m_Yaw,   R.m_Yaw),
            xmath::MinAngleDiff(m_Roll,  R.m_Roll)
        };
    }

    //------------------------------------------------------------------------------
    // Difference
    //------------------------------------------------------------------------------
    //
    // Computes the Euclidean difference (magnitude of vector difference).
    //
    // Parameters:
    //  R - The other radian3.
    //
    // Returns:
    //  The difference as a single radian.
    //
    constexpr radian radian3::Difference(const radian3& R) const noexcept
    {
        return radian{ xmath::Sqrt(xmath::Sqr(m_Pitch.m_Value - R.m_Pitch.m_Value) + Sqr(m_Yaw.m_Value - R.m_Yaw.m_Value) + Sqr(m_Roll.m_Value - R.m_Roll.m_Value)) };
    }

    //------------------------------------------------------------------------------
    // isInrange
    //------------------------------------------------------------------------------
    //
    // Checks if all angles are within the given min and max.
    //
    // Parameters:
    //  Min - The minimum radian.
    //  Max - The maximum radian.
    //
    // Returns:
    //  True if all angles are in range, false otherwise.
    //
    constexpr bool radian3::isInrange(radian Min, radian Max) const noexcept
    {
        return isInRange(m_Pitch, Min, Max) && isInRange(m_Yaw, Min, Max) && isInRange(m_Roll, Min, Max);
    }

    //------------------------------------------------------------------------------
    // isValid
    //------------------------------------------------------------------------------
    //
    // Checks if all angles are valid (not NaN or Inf).
    //
    // Returns:
    //  True if valid, false otherwise.
    //
    inline bool radian3::isValid(void) const noexcept
    {
        return xmath::isValid(m_Pitch.m_Value) && xmath::isValid(m_Yaw.m_Value) && xmath::isValid(m_Roll.m_Value);
    }

    //------------------------------------------------------------------------------
    // hasGimbalLock
    //------------------------------------------------------------------------------
    //
    // Checks if the rotation is in gimbal lock (cos(pitch) approx 0).
    //
    // Returns:
    //  True if in gimbal lock, false otherwise.
    //
    inline bool radian3::hasGimbalLock(void) const noexcept
    {
        assert(isValid() && "Invalid angles in hasGimbalLock");
        return xmath::Abs(xmath::Cos(m_Pitch)) <= flt_tol_v;
    }

    //------------------------------------------------------------------------------
    // Normalize
    //------------------------------------------------------------------------------
    //
    // Normalizes the angles to [-pi, pi) using ModAngle2.
    //
    // Returns:
    //  Reference to self.
    //
    inline radian3& radian3::Normalize(void) noexcept
    {
        assert(isValid() && "Invalid angles before Normalize");
        return ModAngle2();
    }

    //------------------------------------------------------------------------------
    // NormalizeCopy
    //------------------------------------------------------------------------------
    //
    // Returns a normalized copy of the angles to [-pi, pi).
    //
    // Returns:
    //  Normalized radian3.
    //
    inline radian3 radian3::NormalizeCopy(void) const noexcept
    {
        assert(isValid() && "Invalid angles in NormalizeCopy");
        radian3 copy = *this;
        return copy.Normalize();
    }

    //------------------------------------------------------------------------------
    // Clamp
    //------------------------------------------------------------------------------
    //
    // Clamps each angle to [Min, Max].
    //
    // Parameters:
    //  Min - Minimum value for all angles.
    //  Max - Maximum value for all angles.
    //
    // Returns:
    //  Reference to self.
    //
    inline radian3& radian3::Clamp(radian Min, radian Max) noexcept
    {
        assert(isValid() && xmath::isValid(Min.m_Value) && xmath::isValid(Max.m_Value) && "Invalid values in Clamp");
        m_Pitch = xmath::Clamp(m_Pitch, Min, Max);
        m_Yaw   = xmath::Clamp(m_Yaw, Min, Max);
        m_Roll  = xmath::Clamp(m_Roll, Min, Max);
        return *this;
    }

    //------------------------------------------------------------------------------
    // ClampCopy
    //------------------------------------------------------------------------------
    //
    // Returns a clamped copy.
    //
    // Parameters:
    //  Min - Minimum value for all angles.
    //  Max - Maximum value for all angles.
    //
    // Returns:
    //  Clamped copy.
    //
    inline radian3 radian3::ClampCopy(radian Min, radian Max) const noexcept
    {
        assert(isValid() && xmath::isValid(Min.m_Value) && xmath::isValid(Max.m_Value) && "Invalid values in ClampCopy");
        radian3 copy = *this;
        return copy.Clamp(Min, Max);
    }

    //------------------------------------------------------------------------------
    // Lerp
    //------------------------------------------------------------------------------
    //
    // Linearly interpolates between this and other, with angle wrapping.
    //
    // Parameters:
    //  t - Interpolation factor (0 to 1).
    //  Other - The target radian3.
    //
    // Returns:
    //  Interpolated radian3.
    //
    inline radian3 radian3::Lerp(float t, const radian3& Other) const noexcept
    {
        assert(isValid() && Other.isValid() && xmath::isValid(t) && "Invalid values in Lerp");
        return radian3{
            LerpAngle(t, m_Pitch, Other.m_Pitch),
            LerpAngle(t, m_Yaw, Other.m_Yaw),
            LerpAngle(t, m_Roll, Other.m_Roll)
        };
    }

    //------------------------------------------------------------------------------
    // operator ==
    //------------------------------------------------------------------------------
    //
    // Equality comparison.
    //
    // Parameters:
    //  R - The other radian3.
    //
    // Returns:
    //  True if equal, false otherwise.
    //
    constexpr bool radian3::operator==(const radian3& R) const noexcept
    {
        return m_Pitch == R.m_Pitch && m_Yaw == R.m_Yaw && m_Roll == R.m_Roll;
    }

    //------------------------------------------------------------------------------
    // operator +=
    //------------------------------------------------------------------------------
    //
    // Adds another radian3 component-wise.
    //
    // Parameters:
    //  R - The radian3 to add.
    //
    // Returns:
    //  Reference to self.
    //
    inline const radian3& radian3::operator+=(const radian3& R) noexcept
    {
        assert(isValid() && R.isValid() && "Invalid angles in operator+=");
        m_Pitch += R.m_Pitch;
        m_Yaw += R.m_Yaw;
        m_Roll += R.m_Roll;
        return *this;
    }

    //------------------------------------------------------------------------------
    // operator -=
    //------------------------------------------------------------------------------
    //
    // Subtracts another radian3 component-wise.
    //
    // Parameters:
    //  R - The radian3 to subtract.
    //
    // Returns:
    //  Reference to self.
    //
    inline const radian3& radian3::operator-=(const radian3& R) noexcept
    {
        assert(isValid() && R.isValid() && "Invalid angles in operator-=");
        m_Pitch -= R.m_Pitch;
        m_Yaw -= R.m_Yaw;
        m_Roll -= R.m_Roll;
        return *this;
    }

    //------------------------------------------------------------------------------
    // operator *=
    //------------------------------------------------------------------------------
    //
    // Scales all angles by a scalar.
    //
    // Parameters:
    //  Scalar - The scalar to multiply by.
    //
    // Returns:
    //  Reference to self.
    //
    inline const radian3& radian3::operator*=(float Scalar) noexcept
    {
        assert(isValid() && xmath::isValid(Scalar) && "Invalid values in operator*=");
        m_Pitch *= Scalar;
        m_Yaw *= Scalar;
        m_Roll *= Scalar;
        return *this;
    }

    //------------------------------------------------------------------------------
    // operator /=
    //------------------------------------------------------------------------------
    //
    // Divides all angles by a scalar.
    //
    // Parameters:
    //  Scalar - The scalar to divide by.
    //
    // Returns:
    //  Reference to self.
    //
    inline const radian3& radian3::operator/=(float Scalar) noexcept
    {
        assert(isValid() && xmath::isValid(Scalar) && Scalar != 0.0f && "Invalid or zero scalar in operator/=");
        m_Pitch /= Scalar;
        m_Yaw /= Scalar;
        m_Roll /= Scalar;
        return *this;
    }

    //------------------------------------------------------------------------------
    // operator +
    //------------------------------------------------------------------------------
    //
    // Adds another radian3 component-wise.
    //
    // Parameters:
    //  R - The radian3 to add.
    //
    // Returns:
    //  The sum.
    //
    inline radian3 radian3::operator+(const radian3& R) const noexcept
    {
        assert(isValid() && R.isValid() && "Invalid angles in operator+");
        return radian3{ m_Pitch + R.m_Pitch, m_Yaw + R.m_Yaw, m_Roll + R.m_Roll };
    }

    //------------------------------------------------------------------------------
    // operator -
    //------------------------------------------------------------------------------
    //
    // Subtracts another radian3 component-wise.
    //
    // Parameters:
    //  R - The radian3 to subtract.
    //
    // Returns:
    //  The difference.
    //
    inline radian3 radian3::operator-(const radian3& R) const noexcept
    {
        assert(isValid() && R.isValid() && "Invalid angles in operator-");
        return radian3{ m_Pitch - R.m_Pitch, m_Yaw - R.m_Yaw, m_Roll - R.m_Roll };
    }

    //------------------------------------------------------------------------------
    // operator -
    //------------------------------------------------------------------------------
    //
    // Negates the radian3.
    //
    // Returns:
    //  The negated radian3.
    //
    inline radian3 radian3::operator-(void) const noexcept
    {
        assert(isValid() && "Invalid angles in unary operator-");
        return radian3{ -m_Pitch, -m_Yaw, -m_Roll };
    }

    //------------------------------------------------------------------------------
    // operator /
    //------------------------------------------------------------------------------
    //
    // Divides by a scalar.
    //
    // Parameters:
    //  S - The scalar.
    //
    // Returns:
    //  The result.
    //
    inline radian3 radian3::operator/(float S) const noexcept
    {
        assert(isValid() && xmath::isValid(S) && S != 0.0f && "Invalid or zero scalar in operator/");
        return radian3{ m_Pitch / S, m_Yaw / S, m_Roll / S };
    }

    //------------------------------------------------------------------------------
    // operator *
    //------------------------------------------------------------------------------
    //
    // Multiplies by a scalar.
    //
    // Parameters:
    //  S - The scalar.
    //
    // Returns:
    //  The result.
    //
    inline radian3 radian3::operator*(float S) const noexcept
    {
        assert(isValid() && xmath::isValid(S) && "Invalid values in operator*");
        return radian3{ m_Pitch * S, m_Yaw * S, m_Roll * S };
    }

    //------------------------------------------------------------------------------
    // operator *
    //------------------------------------------------------------------------------
    //
    // Multiplies scalar by radian3 (commutative).
    //
    // Parameters:
    //  S - The scalar.
    //  R - The radian3.
    //
    // Returns:
    //  The result radian3.
    //
    inline radian3 operator*(float S, const radian3& R) noexcept
    {
        assert(R.isValid() && isValid(S) && "Invalid values in operator*");
        return R * S;
    }
}