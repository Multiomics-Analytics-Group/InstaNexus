import sys
import argparse
from . import script_dbg, script_greedy

def main():
    banner = r"""
     ______                   __             __  __                                   
    /\__  _\                 /\ \__         /\ \/\ \                                  
    \/_/\ \/     ___     ____\ \ ,_\    __  \ \ `\\ \     __   __  _  __  __    ____  
       \ \ \   /' _ `\  /',__\\ \ \/  /'__`\ \ \ , ` \  /'__`\/\ \/'\/\ \/\ \  /',__\ 
        \_\ \__/\ \/\ \/\__, `\\ \ \_/\ \L\.\_\ \ \`\ \/\  __/\/>  </\ \ \_\ \/\__, `\
        /\_____\ \_\ \_\/\____/ \ \__\ \__/.\_\\ \_\ \_\ \____\/\_/\_\\ \____/\/\____/
        \/_____/\/_/\/_/\/___/   \/__/\/__/\/_/ \/_/\/_/\/____/\//\/_/ \/___/  \/___/ 
    """

    parser = argparse.ArgumentParser(
        prog="instanexus",
        description=(banner + "\n"
            "InstaNexus CLI: de novo protein sequencing based on InstaNovo,\n\n" \
            "an end-to-end workflow from de novo peptides to proteins\n\n"
            "Usage:\n"
            "  instanexus <command> [options]\n\n"
            "Available commands:\n"
            "  dbg       Run De Bruijn Graph assembly pipeline\n"
            "  greedy    Run greedy assembly pipeline\n\n"
            "Examples:\n"
            "  instanexus dbg --input_csv inputs/sample.csv --chain light --folder_outputs outputs --reference\n"
            "  instanexus greedy --input_csv inputs/sample.csv --folder_outputs outputs\n\n"
            "Use 'instanexus <command> --help' for detailed options."
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    
    parser.add_argument('--version', action='version', version='InstaNexus 0.1.0'),

    subparsers = parser.add_subparsers(dest="command", help="subcommands")

    # subcommands
    subparsers.add_parser("dbg", help="Run de Bruijn graph assembly pipeline")
    subparsers.add_parser("greedy", help="Run greedy assembly pipeline")

    args, extra = parser.parse_known_args()

    if args.command == "dbg":
        sys.argv = [sys.argv[0]] + extra
        script_dbg.cli()
    elif args.command == "greedy":
        sys.argv = [sys.argv[0]] + extra
        script_greedy.cli()
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
