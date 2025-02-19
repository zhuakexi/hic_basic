import argparse
# To register a new command, create a new module in the same directory and import it here (1/2).
from . import (
    version,
    render,
)

# To register a new command, append the module to the list below (2/2).
modules = [
    version,
    render,
]
def register_commands(subparsers, modules):
    """
    Register subcommands.
    """
    for module in modules:
        subparser = subparsers.add_parser(module.name, help=module.description)
        module.add_arguments(subparser)
        subparser.set_defaults(func=module.run)

def cli():
    """
    Setup argument parser and dispatch commands to their handlers.
    """
    parser = argparse.ArgumentParser(description="CLI tool for interacting with my project.")
    subparsers = parser.add_subparsers(dest='command', help='Available commands')

    register_commands(subparsers, modules)
    
    args = parser.parse_args()
    
    # If no command is provided, print help and exit
    if args.command is None or args.command.lower() in ['help', '-h', '--help']:
        parser.print_help()
        return
    
    # Ensure that the command is properly registered
    if not hasattr(args, "func"):
        print("Error: Command not recognized or not properly registered.")
        parser.print_help()
        return
    
    # Dispatch the command to its handler
    args.func(args)
    
