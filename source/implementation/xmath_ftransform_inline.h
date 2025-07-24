namespace xmath
{
    inline transform3 transform3::Blend(const transform3& From, const transform3& To, const float T) noexcept
    {
        transform3 B;
        B.m_Scale       = From.m_Scale + T * (To.m_Scale - From.m_Scale);
        B.m_Rotation    = From.m_Rotation.Slerp(To.m_Rotation, T);
        B.m_Position    = From.m_Position + T * (To.m_Position - From.m_Position);
        return B;
    }

    inline transform3& transform3::Blend(const transform3& To, const float T) noexcept
    {
        m_Scale        += T * (To.m_Scale - m_Scale);
        m_Rotation      = m_Rotation.Slerp(To.m_Rotation, T);
        m_Position     += T * (To.m_Position - m_Position);
        return *this;
    }

    inline void transform3::setupIdentity(void) noexcept
    {
        m_Position      = fvec3::fromZero();
        m_Scale         = fvec3::fromOne();
        m_Rotation      = fquat::fromIdentity();
    }

    inline fmat4 transform3::toMatrix(void) const noexcept
    {
        return { m_Position, m_Rotation, m_Scale };
    }

    inline transform3::operator fmat4 () const noexcept
    {
        return { m_Position, m_Rotation, m_Scale };
    }


    inline transform2 transform2::Blend(const transform2& From, const transform2& To, const float T) noexcept
    {
        transform2 B;
        B.m_Scale       = From.m_Scale    + T * (To.m_Scale    - From.m_Scale);
        B.m_Rotation    = From.m_Rotation + T * (To.m_Rotation - From.m_Rotation);
        B.m_Position    = From.m_Position + T * (To.m_Position - From.m_Position);
        return B;
    }

    inline transform2& transform2::Blend(const transform2& To, const float T) noexcept
    {
        m_Scale     += T * (To.m_Scale    - m_Scale);
        m_Rotation  += T * (To.m_Rotation - m_Rotation);
        m_Position  += T * (To.m_Position - m_Position);
        return *this;
    }

    inline void transform2::setupIdentity(void) noexcept
    {
        m_Position  = fvec2::fromZero();
        m_Scale     = fvec2::fromOne();
        m_Rotation  = 0_xdeg;
    }


}