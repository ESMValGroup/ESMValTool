
class Error(Exception):
    """Base class for exceptions in this module."""
    pass


class GeneralError(Error):
    """Exception raised for general errors.

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self, expression, message):
        self.expression = expression
        self.message = message
        print(("GeneralError: " + expression + " / " + message))


class PathError(Error):
    """Exception raised for errors in the input.

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self, expression, message):
        self.expression = expression
        self.message = message
        print(("PathError: " + expression + " / " + message))


class ConfigurationError(Error):
    """Exception raised for errors in the configuration.

    Attributes:
        expression  -- input expression in which the error occurred
        message     -- explanation of the error
    """

    def __init__(self, expression, message):
        self.expression = expression
        self.message = message
        print(("ConfigError: " + expression + " / " + message))


class ImplementationError(Error):
    """ Exception raised for features not implemented yet

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self, expression, message):
        self.expression = expression
        self.message = message
        print(("ImplementationError: " + expression + " / " + message))
        
        
class EmptyContentError(Error): 
    """ Exception raised for empty content

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self, expression, message):
        self.expression = expression
        self.message = message
        print(("contentError: " + expression + " / " + message))


