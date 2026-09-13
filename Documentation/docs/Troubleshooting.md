# Issues and Troubleshooting

The following page includes some of the issues I have faced in using the ECCP program and some of the tips and troubleshooting advice for overcoming various difficulties. 

## ``ImportError`` or ``ModuleNotFoundError`` when working in the folder that holds your clones

If you keep your clones of these programs together in one folder, like this:

~~~
Solar_Cell_Project/
|-- ACSD/
|-- ECCP/
|-- SUMELF/
|-- ...
~~~

then running a command or a python script **from that folder** fails with an error such as:

~~~
ImportError: cannot import name '__version__' from 'ECCP' (unknown location)
~~~

or:

~~~
ModuleNotFoundError: No module named 'SUMELF'
~~~

even though you installed everything correctly.

This happens because python searches the current directory first. The ``ECCP`` folder sitting there is the **git repository**, not the python package: the package is one level further in, at ``ECCP/ECCP``. The repository folder has no ``__init__.py``, so python treats it as an empty namespace package, finds no code inside it, and stops looking. Your installed copy is never reached.

The fix is to work from anywhere other than that folder. Change into the directory holding the crystals you are working on and run from there, which is what you would normally be doing anyway:

~~~bash
cd /path/to/my_crystals
eccp --help
~~~

!!! tip

	Only the folder that *directly* contains the repository folders is affected. Subfolders of it are fine, and so is any unrelated directory.

## Other Issues

This program is definitely a "work in progress". I have made it as easy to use as possible, but there are always oversights to program development and some parts of it may not be as easy to use as it could be. 

If you have any issues with the program or you think there would be better/easier ways to use and implement things in the ECCP program, write a message on the [ECCP Github Issues page](https://github.com/geoffreyweal/ECCP/issues). Feedback is very much welcome!
