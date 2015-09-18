if __name__ == "phe":
    # If this package is added as library, append extended path.
    from pkgutil import extend_path
    __path__ = extend_path(__path__, __name__)