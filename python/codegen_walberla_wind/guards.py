def require_found(found, description, package):
    if found is None or found is False:
        raise ValueError(f"{package} no longer provides {description}")
    return found
