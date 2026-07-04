#pragma once

namespace fred {

// Abstract interface for coolant thermophysical properties used by subchannel models.
class SubchannelCoolantProperties {
public:
    virtual double cp (double T) const = 0;   // specific heat     [J/(kg·K)]
    virtual double rho(double T) const = 0;   // density           [kg/m³]
    virtual double k  (double T) const = 0;   // conductivity      [W/(m·K)]
    virtual double mu (double T) const = 0;   // dynamic viscosity [Pa·s]
    virtual ~SubchannelCoolantProperties() = default;
};

} // namespace fred
