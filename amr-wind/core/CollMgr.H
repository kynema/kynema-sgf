#ifndef COLLMGR_H
#define COLLMGR_H

#include <string>
#include <memory>
#include <unordered_map>
#include "AMReX_Vector.H"

namespace amr_wind {

/** Interface to manage objects of sub-classes based on factory registration
 *  \ingroup core
 *
 *  Provides boiler-plate code to register a collection and then query using
 *  lookup keywords
 *
 *  \tparam Collection Class holding instances of collection (e.g., PhysicsMgr)
 *  \tparam Base class of the collection (e.g., Physics)
 */
template <typename Collection, typename Type, bool allow_duplicates = false>
class CollMgr
{
public:
    using TypePtr = std::unique_ptr<Type>;
    using TypeVector = amrex::Vector<TypePtr>;

private:
    CollMgr() = default;

public:
    ~CollMgr() = default;

    CollMgr(const CollMgr&) = delete;
    CollMgr& operator=(const CollMgr&) = delete;

    /** Create a new object and register it in the vector and setup lookup
     * options
     */
    template <class... Args>
    Type& create(const std::string& key, Args&&... args)
    {
        if (!allow_duplicates && contains(key)) {
            return operator()(key);
        }

        m_obj_vec.emplace_back(Type::create(key, std::forward<Args>(args)...));
        m_obj_map[key] = m_obj_vec.size() - 1;

        return *m_obj_vec.back();
    }

    //! Return a vector of the registered objects
    TypeVector& objects() { return m_obj_vec; }
    [[nodiscard]] const TypeVector& objects() const { return m_obj_vec; }

    //! Query if an object exists using the lookup key
    [[nodiscard]] bool contains(const std::string& key) const
    {
        auto it = m_obj_map.find(key);
        return (it != m_obj_map.end());
    }

    //! Return the object corresponding to a lookup key
    Type& operator()(const std::string& key)
    {
        return *m_obj_vec[m_obj_map.at(key)];
    }

    const Type& operator()(const std::string& key) const
    {
        return *m_obj_vec[m_obj_map.at(key)];
    }

    //! Return object for a lookup key cast into its exact class definition
    template <typename T>
    T& operator()(const std::string& key)
    {
        return dynamic_cast<T&>(operator()(key));
    }

protected:
    //! Collection of objects registered so far
    TypeVector m_obj_vec;

    //! Key word based lookup
    std::unordered_map<std::string, int> m_obj_map;
    friend Collection;
};

} // namespace amr_wind

#endif /* COLLMGR_H */
