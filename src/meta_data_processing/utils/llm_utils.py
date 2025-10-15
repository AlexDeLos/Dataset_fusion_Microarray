from langchain_core.prompts import ChatPromptTemplate
from langchain_core.messages import SystemMessage
from langchain_core.output_parsers import PydanticOutputParser, StrOutputParser
from langchain.chat_models import init_chat_model
from typing import List, Literal
from pydantic import ValidationError, BaseModel, Field


class MetaData(BaseModel):
    tissue: Literal["root", "leaf", "flower", "shoot", "rosette", "bud", "whole_plant", "silique", "callus", "seed", "seedling", "unknown"] = Field(..., description="Tissue the samples was extracted from.") # "stem", 
    # tissue: str = Field(..., description="Tissue the samples was extracted from.") # "stem", 

    # Field with description explaining its purpose
    treatment: List[Literal[
    "Drought Stress",
    "Salinity Stress",
    "Heat Stress",
    "Cold Stress",
    "Chemical Stress",
    "Nutrient Deficiency",
    "Pathogen Attack",
    "Low Light Stress",
    "High Light Stress",
    "Red Light Stress",
    "Other Light Stress",
    "Other stress",
    "No stress"]] = Field(..., description="List of treatments and stresses that was applied to the sample, each unique stress or treatment should have one, and only one, entry in this list")
    
    # Field with description explaining its purpose
    medium: str = Field(..., description="Growth medium of the sample.")
    # Field with description explaining its purpose
    # age: float = Field(..., description="The age of the plant since germination in days.") 
    # # Field with description explaining its purpose
    # treatment_time: float = Field(..., description="Time in hours that between the application of the treatment to the sample being harvested.") #! this is not as accurate




def get_metadata(study_info:dict,sample_info:dict,model:str='gemini-2.5-flash',temp:float=0):
    llm = init_chat_model(model=model,
                        model_provider="google_genai",
                        temperature=temp)
    parser = PydanticOutputParser(pydantic_object=MetaData)
    system_message = SystemMessage(content="You are a biology reseracher that has been tasked with creating a script for analysing the metadata of given samples that follows a particular sctructure and returning information on the biological experimental conditions of the sample that are requested. You want to stick to labels that are both relevant but wide enouogh that other studies can fall under the same label. If more than one stress is applied return both.")


    format_instructions = parser.get_format_instructions()

    prompt = ChatPromptTemplate.from_messages([
        ("system", "Extract per schema:\n{format_instructions}"),
        ("human", "{text}"),
    ]).partial(format_instructions=format_instructions)

    parsing_llm = prompt | llm | parser


    parsing_prompt = '''
    <task>
    Extract the relevant biological experimental information of the sample following the given schema:
    </task>
    <metadata>
        <study_metadata>
            {study_info}
        </study_metadata>
        <sample_metadata>
            {sample_info}
        </sample_metadata>
    </metadata>
    '''    
    try:
        result = parsing_llm.invoke({"text": parsing_prompt.format(study_info=study_info,sample_info=sample_info)}) # type: ignore #! may also raise an uncaght JSONDecodeError ->GSM1126269
    except Exception as e:
        print(f'trying a better model-> {e}')
        if model=='gemini-2.5-pro':
            raise e#-> GSM1126262 ->GSM843683
        return get_metadata(study_info,sample_info, model='gemini-2.5-pro', temp=0.1)
    return result

def get_metadata_script(study_info:dict,sample_info:dict, study_id: str,model:str='gemini-2.5-flash',temp:float=0):
    llm = init_chat_model(model=model,
                        model_provider="google_genai",
                        temperature=temp)
    parser = PydanticOutputParser(pydantic_object=MetaData)


    format_instructions = parser.get_format_instructions()

    prompt = ChatPromptTemplate.from_messages([
        ("system", "Produce a python cade function that can extract the following schema from the two data dictionary objects:\n{format_instructions}"),
        ("human", "{text}"),
    ]).partial(format_instructions=format_instructions)

    parsing_llm = prompt | llm | StrOutputParser()


    parsing_prompt = '''
    <task>
    Provide a python function function with signature (def {study_id}_extractor(sample_metadata: dict) -> dict) that extracts the relevant biological experimental information of the sample from this dataset following the given schema.
    When the code is ran it should give proper labels for (following the schema) for the sample that is represented in the sample_metadata.
    </task>
    <guidance>
        Keep in mind that information like meduim or tissue can be constant across samples.
        Focus on location WHERE in the dictionary the desired information is stored.
    </guidance>
    <metadata>
        The following is the study metadata (constant across all samples in the study) and list of the metadata of the samples in this studies.
        Use it to guide to know where in the sample_metadata to find the information.
        <study_metadata>
            {study_info}
        </study_metadata>
        <samples_metadata>
            {sample_info}
        </samples_metadata>
    </metadata>
    '''

    try:
        result = parsing_llm.invoke({"text": parsing_prompt.format(study_id=study_id,study_info=study_info,sample_info=sample_info)}) # type: ignore #! may also raise an uncaght JSONDecodeError ->GSM1126269
    except Exception as e:
        print(f'trying a better model-> {e}')
        if model=='gemini-2.5-pro':
            raise e#-> GSM1126262 ->GSM843683
        result =  get_metadata_script(study_info,sample_info, study_id,'gemini-2.5-pro')
    return result




def get_condensed_labels(study_info:dict,sample_info:dict,model:str='gemini-2.5-flash',temp:float=0):
    llm = init_chat_model(model=model,
                        model_provider="google_genai",
                        temperature=temp)
    parser = PydanticOutputParser(pydantic_object=MetaData)


    format_instructions = parser.get_format_instructions()

    prompt = ChatPromptTemplate.from_messages([
        ("system", "You are a biology reseracher that has been tasked with creating a script for analysing the metadata of given samples that follows a particular sctructure and returning information on the biological experimental conditions of the sample that are requested. You want to stick to labels that are both relevant but wide enouogh that other studies can fall under the same label. If more than one stress is applied label it as 'control'. Keep the labels in the following schema:\n{format_instructions}"),
        ("human", "{text}"),
    ]).partial(format_instructions=format_instructions)

    parsing_llm = prompt | llm | parser


    parsing_prompt = '''
    <task>
    Using the study information and the extracted sample labels ground the labels to comonly used term matching the ontology:
    </task>
    <metadata>
        <study_metadata>
            {study_info}
        </study_metadata>
        <samples_label>
            {sample_info}
        </samples_label>
    </metadata>
    '''    
    try:
        result = parsing_llm.invoke({"text": parsing_prompt.format(study_info=study_info,sample_info=sample_info)}) # type: ignore #! may also raise an uncaght JSONDecodeError ->GSM1126269
    except Exception as e:
        print(f'trying a better model-> {e}')
        if model=='gemini-2.5-pro':
            raise e#-> GSM1126262 ->GSM843683
        return get_metadata(study_info,sample_info, model='gemini-2.5-pro', temp=0.1)
    return result