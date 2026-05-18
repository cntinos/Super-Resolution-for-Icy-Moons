<style>body {text-align: justify}</style>
# Super-Resolution for Icy Moons

## Description
The icy moons of the outer solar system have been found in the forefront of planetary exploration, due to their diverse geomorphology coupled with numerous discoveries hinting to possible habitable conditions. Of paramount importance for constraining the internal structure and properties of these bodies are high-precision cartographic products, allowing for linking surface morphologies and composition to internal activity. These cartographic products are constructed by reconciling numerous different observations, which, due to the complexities of an orbital tour, are recorded under a great range of different illumination and observation conditions, with varying spatial sampling. The aim of this study is to provide a complete pipeline for producing Super-Resolution (SR) maps for icy moons, by combining all the overlapping observations for a specific target. The development and testing of our pipeline are carried out for Enceladus, using the hyperspectral cubes recorded by the Visual and Infrared Mapping Spectrometer (VIMS) instrument, on board the Cassini Saturn orbiter. We first navigate the data, calculating accurate projections that can be combined into global mosaics of Enceladus. We then implement a global photometric correction model that compensates for the different acquisition conditions. Finally, we combine all the overlapping observations to calculate the new SR maps of Enceladus in the infrared. The new SR global map of Enceladus provides new insight into the moon’s geomorphology as revealed from its surface spectral heterogeneity at finer spatial resolution, demonstrating the capability of SR to augment our current understanding of icy moons and help prepare for future observations.

## Badges
On some READMEs, you may see small images that convey metadata, such as whether or not all the tests are passing for the project. You can use Shields to add some to your README. Many services also have instructions for adding a badge.

## Pipeline

![Pipeline flowchart](sr_pipeline_flowchart.svg)

## Installation
Within a particular ecosystem, there may be a common way of installing things, such as using Yarn, NuGet, or Homebrew. However, consider the possibility that whoever is reading your README is a novice and would like more guidance. Listing specific steps helps remove ambiguity and gets people to using your project as quickly as possible. If it only runs in a specific context like a particular programming language version or operating system or has dependencies that have to be installed manually, also add a Requirements subsection.

## Usage
Use examples liberally, and show the expected output if you can. It's helpful to have inline the smallest example of usage that you can demonstrate, while providing links to more sophisticated examples if they are too long to reasonably include in the README.

## Support
Tell people where they can go to for help. It can be any combination of an issue tracker, a chat room, an email address, etc.

## Roadmap
If you have ideas for releases in the future, it is a good idea to list them in the README.

## Contributing
State if you are open to contributions and what your requirements are for accepting them.

For people who want to make changes to your project, it's helpful to have some documentation on how to get started. Perhaps there is a script that they should run or some environment variables that they need to set. Make these steps explicit. These instructions could also be useful to your future self.

You can also document commands to lint the code or run tests. These steps help to ensure high code quality and reduce the likelihood that the changes inadvertently break something. Having instructions for running tests is especially helpful if it requires external setup, such as starting a Selenium server for testing in a browser.

## Authors and acknowledgment
Show your appreciation to those who have contributed to the project.

## License
For open source projects, say how it is licensed.

## Project status
If you have run out of energy or time for your project, put a note at the top of the README saying that development has slowed down or stopped completely. Someone may choose to fork your project or volunteer to step in as a maintainer or owner, allowing your project to keep going. You can also make an explicit request for maintainers.
