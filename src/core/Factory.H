#ifndef FACTORY_H
#define FACTORY_H

#include <memory>
#include <unordered_map>
#include <iostream>

#include "AMReX_Print.H"
#include "AMReX.H"

namespace amr_wind {

/** Runtime selection of options based on string keyword lookup
 *  \ingroup core
 *
 *  Factory provides an automated way to register subclasses that can be
 *  instantiated during runtime using string based lookups. The facility mimics
 *  the runTimeSelectionTable in OpenFOAM and the current structure is adapted
 *  from http://www.nirfriedman.com/2018/04/29/unforgettable-factory/
 *
 *  The usage is described using the following example from code. We would like
 *  to support creation of a turbulence model based on user inputs that are read
 *  from a file. We will create a base class TurbulenceModel that defines the
 *  API for the turbulence model and how it interacts with the rest of the code.
 *  Then new models, e.g., Smagorinsky, k-Omega SST/SAS, etc., are added as
 *  subclasses to this instance.
 *
 *  ```
 *  // Create a base class that can be used to register and create instances
 *  //
 *  // The template arguments are the classname for CRTP and types of additional
 *  // arguments to be passed to the constructor of the object
 *  class TurbulenceModel : public Factory<TurbulenceModel, const CFDSim&>
 *  {
 *       // Define a base identifier string for debugging purposes
 *       static const std::string base_identifier() { return "TurbulenceModel";
 * }
 *
 *       // Define API as appropriate here
 *  };
 *
 *  // Define a subclass
 *  // Instead of using the base class we use Base::Register to allow
 * registration class Smagorinsky : public
 * TurbulenceModel::Register<Smagorinsky>
 *  {
 *       // This is the keyword that is used to lookup and create this instance
 *       static const std::string identifier() { return "Smagorinsky"; }
 *
 *       // Model implementation here
 *  };
 *
 *  // To create a turbulence model instance
 *  auto sim = CFDSim(mesh);
 *  // tmodel is a std::unique_ptr<TurbulenceModel> instance
 *  auto tmodel = TurbulenceModel::create("Smagorinsky", sim);
 *  ```
 *
 *  \tparam Base base class (e.g., Physics)
 *  \tparam Args Arguments types to be passed to the constructor during
 * instantiation
 */
template <class Base, class... Args>
struct Factory
{
    /** Create an instance of the concrete subclass based on the runtime keyword
     *
     *  \param key Unique identifier used to lookup subclass instance
     *  \param args Arguments to the constructor of the subclass
     */
    static std::unique_ptr<Base> create(const std::string& key, Args... args)
    {
        key_exists_or_error(key);
        auto ptr = table().at(key)(std::forward<Args>(args)...);
        amrex::Print() << "Creating " << Base::base_identifier()
                       << " instance: " << key << std::endl;
        return ptr;
    }

    /** Print a list of the valid subclasses registered to this base instance
     *
     *  \param os Valid output stream
     */
    static void print(std::ostream& os)
    {
        const auto& tbl = table();
        os << Base::base_identifier() << " " << tbl.size() << std::endl;
        for (const auto& it : tbl) {
            os << " - " << it.first << std::endl;
        }
    }

    //! Class to handle registration of subclass for runtime selection
    template <class T>
    struct Register : public Base
    {
        friend T;
        static bool registered;

        static bool add_sub_type()
        {
            // TODO: Add checks against multiple registration
            Factory::table()[T::identifier()] =
                [](Args... args) -> std::unique_ptr<Base> {
                return std::unique_ptr<Base>(
                    new T(std::forward<Args>(args)...));
            };
            return true;
        }

        ~Register() override
        {
            if (registered) {
                const auto& tbl = Factory::table();
                const auto it = tbl.find(T::identifier());
                registered = (it != tbl.end());
            }
        }

    private:
        Register() { (void)registered; }
    };

    virtual ~Factory() = default;
    friend Base;

private:
    using CreatorFunc = std::unique_ptr<Base> (*)(Args...);
    using LookupTable = std::unordered_map<std::string, CreatorFunc>;

    //! Check if the keyword exists in table or print diagnostic message and
    //! abort
    static void key_exists_or_error(const std::string& key)
    {
        const auto& tbl = table();
        if (tbl.find(key) == tbl.end()) {
            // Print valid options
            if (amrex::ParallelDescriptor::IOProcessor()) {
                print(std::cout);
            }
            // Quit with error
            amrex::Abort(
                "In " + Base::base_identifier() +
                " cannot find instance: " + key);
        }
    }

    //! Lookup table containing all registered instances
    static LookupTable& table()
    {
        static LookupTable tbl;
        return tbl;
    }

    Factory() = default;
};

template <class Base, class... Args>
template <class T>
bool Factory<Base, Args...>::Register<T>::registered =
    Factory<Base, Args...>::Register<T>::add_sub_type();

} // namespace amr_wind

#endif /* FACTORY_H */
