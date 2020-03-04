#ifndef PHASETRACER_PROPERTY_HPP_INCLUDED
#define PHASETRACER_PROPERTY_HPP_INCLUDED
/**
   Convenience macro for adding getters and setters to our objects
*/

#define PROPERTY(type, name, default) \
 public: \
  void set_##name(type _##name) {name = _##name;} \
  type get_##name() const {return name;} \
 private: \
  type name = default;

#define PROTECTED_PROPERTY(type, name, default) \
 public: \
  void set_##name(type _##name) {name = _##name;} \
  type get_##name() const {return name;} \
 protected: \
  type name = default;

#define PROPERTY_CUSTOM_SETTER(type, name, default) \
 public: \
  type get_##name() const {return name;} \
 private: \
  type name = default; \

#define PROTECTED_PROPERTY_CUSTOM_SETTER(type, name, default) \
 public: \
  type get_##name() const {return name;} \
 protected: \
  type name = default;

#endif
