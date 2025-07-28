import yaml
from dataclasses import dataclass
from copy import deepcopy


# Register a tuple constructor with PyYAML
def tuple_constructor(loader, node):
    """Load the sequence of values from the YAML node and returns a tuple constructed from the sequence."""
    values = loader.construct_sequence(node)
    return tuple(values)


yaml.SafeLoader.add_constructor('tag:yaml.org,2002:python/tuple', tuple_constructor)


@dataclass
class ParametersBase:
    """Base class for parameter dataclassess
    """
    section_name: str = "default"
    is_space_to_earth: bool = False  # whether the system is a space station or not

    # whether to enable recursive parameters setting on .yaml file
    # TODO: make every system have this be True and remove this attribute
    nested_parameters_enabled: bool = False

    # TODO: make this be directly called as the default load method, after reading .yml file
    def load_subparameters(self, ctx: str, params: dict, quiet=True):
        """
            ctx: provides information on what subattribute is being parsed.
                 This is mainly for debugging/logging/error handling
            params: dict that contains the attributes needed by the
            quiet: if True the parser will not warn about unset paramters
        """
        # Load all the parameters from the configuration file
        attr_list = [
            a for a in dir(self) if not a.startswith('_') and not callable(getattr(self, a)) and a not in
                    ["section_name", "nested_parameters_enabled",]]

        params_keys = params.keys()

        for k in params_keys:
            if k not in attr_list:
                raise ValueError(
                    f"The parameter {ctx}.{k} was passed, but it doesn't exist on parameters definitions!"
                )

        for attr_name in attr_list:
            default_attr_value = getattr(self, attr_name)

            if attr_name not in params:
                if not quiet:
                    print(
                        f"[INFO]: WARNING. Using default parameters for {ctx}.{attr_name}: {default_attr_value}",
                    )
            elif isinstance(default_attr_value, ParametersBase):
                if not isinstance(params[attr_name], dict):
                    raise ValueError(
                        f"ERROR: Cannot parse section {ctx}.{attr_name}, is {params[attr_name]} instead of a dictionary",
                    )

                # try to recursively set config
                # is a bit hacky and limits some stuff, since it doesn't know the context it is in
                # for example, it cannot get system frequency to set some value
                default_attr_value.load_subparameters(
                    f"{ctx}.{attr_name}", params[attr_name],
                )
            elif isinstance(default_attr_value, list):
                if not isinstance(params[attr_name], list):
                    raise ValueError(
                        f"ERROR: Cannot parse parameter {ctx}.{attr_name}, is \
                            {params[attr_name]} instead of a list",
                    )
                loaded_attr_vals = list()
                default_item = default_attr_value[0]

                if isinstance(default_item, ParametersBase):
                    for prm in params[attr_name]:
                        new_item = deepcopy(default_item)
                        new_item.load_subparameters(
                            f"{self.section_name}.{attr_name}", prm,
                        )
                        loaded_attr_vals.append(new_item)
                else:
                    for prm in params[attr_name]:
                        if not isinstance(prm, type(default_item)):
                            raise ValueError(
                                f"ERROR: Cannot parse section {ctx}.{attr_name}\n"
                                f"List item does not respect expected type of {type(default_item)}\n"
                                f"{prm} has type of {type(prm)}"
                            )
                        loaded_attr_vals.append(prm)

                setattr(self, attr_name, loaded_attr_vals)
            else:
                setattr(self, attr_name, params[attr_name])

    def validate(self, ctx: str):
        """
            This method exists because there was a need to separate params parsing from validating,
                since nested parameters may need some attributes to be set by a "parent"
                before proper validation.
            ctx: context string. It should be a string that gives more information about where
                validation is being called on, so that errors may be better handled/alerted
        """
        attr_list = [
            a for a in dir(self) if not a.startswith('_') and isinstance(getattr(self, a), ParametersBase)
        ]

        for attr in attr_list:
            getattr(self, attr).validate(f"{ctx}.{attr}")

    def load_parameters_from_file(self, config_file: str, quiet=True):
        """Load the parameters from file.
        The sanity check is child class reponsibility

        Parameters
        ----------
        file_name : str
            the path to the configuration file
        quiet: if True the parser will not warn about unset paramters

        Raises
        ------
        ValueError
            if a parameter is not valid
        """

        with open(config_file, 'r') as file:
            config = yaml.safe_load(file)

        if self.section_name.lower() not in config.keys():
            if not quiet:
                print(f"ParameterBase: section {self.section_name} not in parameter file.\
                    Only default parameters where loaded.")
            return

        # Load all the parameters from the configuration file
        attr_list = [
            a for a in dir(self) if not a.startswith('_') and not
            callable(getattr(self, a)) and a != "section_name"
        ]

        params_keys = config[self.section_name].keys()

        for k in params_keys:
            if k not in attr_list:
                raise ValueError(
                    f"The parameter {self.section_name}.{k} was passed, but it doesn't exist on parameters definitions!"
                )

        for attr_name in attr_list:
            try:
                attr_val = getattr(self, attr_name)
                # since there are no tuple types in yaml:
                if isinstance(attr_val, tuple):
                    # TODO: test this conditional. There are no test cases for tuple, and no tuple in any of current parameters
                    # Check if the string defines a list of floats
                    try:
                        param_val = config[self.section_name][attr_name]
                        tmp_val = list(map(float, param_val.split(",")))
                        setattr(self, attr_name, tuple(tmp_val))
                    except ValueError:
                        # its a regular string. Let the specific class implementation
                        # do the sanity check
                        print(
                            f"ParametersBase: could not convert string to tuple \"{self.section_name}.{attr_name}\"",
                        )
                        exit()

                # TODO: make every parameters use this way of setting its own attributes, and remove
                # attr_val.nested_parameters_enabled we check for attr_val.nested_parameters_enabled because we don't
                # want to print notice for this kind of parameter YET
                elif isinstance(attr_val, ParametersBase):
                    if not self.nested_parameters_enabled:
                        continue

                    if not isinstance(config[self.section_name][attr_name], dict):
                        raise ValueError(
                            f"ERROR: Cannot parse section {self.section_name}.{attr_name}, is \
                                {config[self.section_name][attr_name]} instead of a dictionary",
                        )

                    # try to recursively set config
                    # is a bit hacky and limits some stuff, since it doesn't know the context it is in
                    # for example, it cannot get system frequency to set some value
                    attr_val.load_subparameters(
                        f"{self.section_name}.{attr_name}", config[self.section_name][attr_name],
                    )
                elif isinstance(attr_val, list):
                    if not self.nested_parameters_enabled:
                        continue
                    if not isinstance(config[self.section_name][attr_name], list):
                        raise ValueError(
                            f"ERROR: Cannot parse section {self.section_name}.{attr_name}, is \
                                {config[self.section_name][attr_name]} instead of a list",
                        )
                    loaded_attr_vals = list()
                    default_item = attr_val[0]
                    for params in config[self.section_name][attr_name]:
                        new_item = deepcopy(default_item)
                        new_item.load_subparameters(
                            f"{self.section_name}.{attr_name}", params,
                        )
                        loaded_attr_vals.append(new_item)
                    setattr(self, attr_name, loaded_attr_vals)

                else:
                    setattr(
                        self, attr_name,
                        config[self.section_name][attr_name],
                    )

            except KeyError:
                if not quiet:
                    print(
                        f"ParametersBase: NOTICE! Configuration parameter \"{self.section_name}.{attr_name}\" \
                            is not set in configuration file. Using default value {attr_val}",
                    )
